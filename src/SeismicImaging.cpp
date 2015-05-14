/*  =================================================================================
  *  Name        : SeismicImaging.cpp
 *  Author      : Mafijul Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Execute methods declared in SeismicImaging.hpp file
 *  Date        : Feb 19, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#include <cmath>
#include <algorithm>
using::std::sort;
#include <stdlib.h>

#include <iostream>
using std::cout;
using std::ofstream;
using std::endl;
#include <fstream>
#include <string>

#include <complex.h>
#include <fftw3.h>

#include "DataModelParams.hpp"
#include "FreqWavenumberParams.hpp"
#include "RefVelParams.hpp"
#include "SeismicImaging.hpp"

#define PI 3.1415926535897932
#define THRESHVEL 1.0       /// Maximum velocity difference to perform interpolation

/* Core method to generate seismic image using PSPI algorithm *****************/
void SeismicImaging::Imaging(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp)
{
    int chkDataFile, chkVelFile;

    /* Seismic data memory allocation *****************************************/
    dmp->data = alloc2DPointer(dmp->nt, dmp->nx, dmp->data);

    /* Velocity model memory allocation ***************************************/
	dmp->vModel = alloc2DPointer(dmp->nz, dmp->nx, dmp->vModel);

	/* Phase shifted F-K data memory allocation *******************************/
	fwp->fzData = alloc2DComPointer(rvp->maxRefVel, dmp->nx, fwp->fzData);

	/* Read a binary file (Seismic data) **************************************/
	chkDataFile = readBinFile(dmp->dataFileName, dmp->nt, dmp->nTraces, dmp->data);

	/* Read binary file (Velocity model) **************************************/
	chkVelFile = readBinFile(dmp->velFileName, dmp->nz, dmp->nx, dmp->vModel);

	for(int i=0;i<dmp->nz; i++)
		for(int j=0;j<dmp->nx; j++)
			dmp->vModel[i][j] = dmp->vModel[i][j]/2;

	if(chkDataFile == 1 && chkVelFile == 1)
	{
		/*Allocation memory for the seismic image  ********************************/
		dmp->image = alloc2DPointer(dmp->nz, dmp->nx, dmp->image);

		/*Image initialisation ************************************************/
		for(int i=0; i<dmp->nz; i++)
			for(int j=0; j<dmp->nx; j++)
				dmp->image[i][j] = 0;

		/*Compute temporal frequency and spatial wave-number axes *************/
		computeFreqWavenumber(fwp, dmp);

		/*Perform FFT in time sampling direction. i.e., p(x,z,t) ----> P(x,z,w)*/
		tXtoFx(dmp, fwp);

		if(dmp->migType == PSPI)
		{
			/*Compute reference velocities for migration method **************/
			computeRefVel(rvp, dmp);

			/*Method to implement downward continuation ***************************/
			downwardExtrapolation(dmp, fwp, rvp);
		}

		if(dmp->migType == GAZDAG)
			/*Method to implement downward continuation ***************************/
			downwardExtrapolation(dmp, fwp, rvp);

		/*Write the image file ************************************************/
		writeBinFile(dmp->imageFileName, dmp->nz, dmp->nx, dmp->image);
	}
	else
		std::cout << "Please check whether the path of binary files are correct or not. " <<std::endl;
}

/*  2D complex float dynamic memory allocation with n1 X n2 size ***************/
inline float complex ** SeismicImaging::alloc2DComPointer(int n1, int n2, float complex **data)
{
	data = (float complex **) fftw_malloc(sizeof(float complex *) * n1);

	for(int i=0; i<n1; i++)
		data[i] = (float complex *) fftw_malloc(sizeof(float complex) * n2);

	return data;
}

/*  2D float dynamic memory allocation with n1 X n2 size *******************************************/
inline float ** SeismicImaging::alloc2DPointer(int n1, int n2, float **data)
{
	data = (float **) malloc(sizeof(float *) * n1);
	for(int i=0; i<n1; i++)
		data[i] = (float *) malloc(sizeof(float) * n2);
	return data;
}

/*Read binary file and store the information into corresponding data pointer ************************/
inline int SeismicImaging::readBinFile(std::string fileName,  int n1, int n2, float **data )
{
    int fsize, chkFile = 0;
    FILE *fp;
    fp = fopen(fileName.c_str(), "rb");
    float *temp;

    if(!fp)
    	std::cout << "File could not opened. " <<std::endl;
    else
    {
    	fsize = n1 * n2;
    	temp = (float *) malloc(sizeof(float) * fsize);
    	fread(temp,sizeof(float), fsize, fp);
    	chkFile = 1;
    	for(int i=0; i<n2; i++)
    		for(int j=0; j<n1; j++)
    			data[j][i] = temp[i*n1+j];

    	fclose(fp);
    	delete [] temp;
    }
    return chkFile;
}

/* Compute reference velocities for PSPI migration method *****************/
inline void SeismicImaging::computeRefVel(RefVelParams* rvp, DataModelParams* dmp)
{
	float *tempVelocity;
	float minVel, maxVel, dv;

	tempVelocity = (float *) malloc(sizeof(float) * (dmp->nx));

	rvp->refVel = alloc2DPointer(dmp->nz, rvp->maxRefVel, rvp->refVel);

	for(int i=0;i<dmp->nz; i++)
	{
		for(int j=0;j<dmp->nx; j++)
			tempVelocity[j] = dmp->vModel[i][j];

		/* To find the mix and max of velocity array ***********************/
		std::sort(tempVelocity, tempVelocity+dmp->nx);
		minVel = tempVelocity[0];
		maxVel = tempVelocity[dmp->nx-1];
		dv = (maxVel - minVel)/(rvp->maxRefVel-1);
		for(int k=0; k<rvp->maxRefVel; k++)
			rvp->refVel[i][k] = minVel + dv*k;
	}
}

/* Compute frequency and spatial wave-number axis  *************************/
inline void SeismicImaging::computeFreqWavenumber(FreqWavenumberParams* fwp, DataModelParams* dmp)
{
	  float dk, dw;
	  /* Determine the temporal frequency axis */
	  dw = 2*PI/(dmp->dt*dmp->nw);
	  fwp->wT = (float *) malloc(sizeof(float) * (dmp->nw/2+1));
	  fwp->wT[0] = 0.0;
	  for(int i=1; i<(dmp->nw/2+1); i++)
		  fwp->wT[i] = dw*i;

	  /* Determine the spatial wave-number axis; Be careful a little bit here. */
	  dk = 2*PI/(dmp->dx*dmp->nk);
	  fwp->wK = (float *) malloc(sizeof(float) * (dmp->nk));
	  fwp->wK[0] = 0.0;
	  for(int i=1; i<dmp->nk; i++)
		fwp->wK[i] = (i<=dmp->nk/2)?(dk)*i:-dk*(dmp->nk-i);
}

/* Perform fft in temporal direction. p(x,z,t) ----> P(x,z,w) ********************* */
inline void SeismicImaging::tXtoFx(DataModelParams* dmp, FreqWavenumberParams* fwp)
{
	  fftw_plan p1;
	  fftw_complex *out1;
	  double *in;

	  out1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (dmp->nw/2+1));
	  in = (double *) fftw_malloc(sizeof(double) * (dmp->nt));
	  fwp->fxData = alloc2DComPointer(dmp->nw/2+1, dmp->nx, fwp->fxData);

	   /**
		** Very important to make a plan to execute. This plan perform FFT
		** of real input data and the output is complex data. The length of
		** output data is len(input)/2+1.
		**/
	   p1 = fftw_plan_dft_r2c_1d(dmp->nw, in, out1, FFTW_ESTIMATE);
	   for(int i=0; i<dmp->nTraces; i++ )
	   {
		   for(int j=0; j<dmp->nt; j++ )
			   in[j] = dmp->data[j][i];

		   fftw_execute(p1);
		   for(int k=0; k<(dmp->nw/2+1); k++)
			   fwp->fxData[k][i] = out1[k];
	   }

	   fftw_destroy_plan(p1);
	   delete [] in;
	   delete [] out1;
}

inline void SeismicImaging::downwardExtrapolation(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp)
{
      int w1, w2, velIndx1, velIndx2;
      float complex shift;
      fftw_complex *out2, *out1, *in1, *in2;
      double dv, w_v2, Kx2, w_v, sumVel, avgVel;
      double eps = 0.1;
      fftw_plan p1, p2;

      w1 = (int)fwp->minF*dmp->dt*dmp->nw +1;
      w2 = (int)fwp->maxF*dmp->dt*dmp->nw +1;

      out1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (dmp->nk));
      in1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (dmp->nk));
      out2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (dmp->nk));
      in2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (dmp->nk));

      p1 = fftw_plan_dft_1d(dmp->nk, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
      p2 = fftw_plan_dft_1d(dmp->nk, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

      /******** Iterate for every frequency component *********************************/
      for(int i=w1; i<=w2; i++)
      {
         printf("Frequency: %f\n", (i-1)/(dmp->dt*dmp->nw));
         for(int k=0;k<dmp->nx; k++)
			in1[k] = fwp->fxData[i][k];

         /****** Iterate for every depth (dz) slice of the velocity model *************/
         for(int j=0; j<dmp->nz;  j++)
         {
        	 fftw_execute(p1);

        	 for(int k=0;k<dmp->nx; k++)
				 /*if(dmp->vModel[j][k] > 0.0)
					 in1[k] = fwp->fxData[i][k] * cexp(dmp->dz * (fwp->wT[i]/dmp->vModel[j][k])*I);
				 else*/
					 in1[k] = fwp->fxData[i][k];

    	     if(dmp->migType == GAZDAG)
    	     {
    	    	 /********* Phase shift due to velocity *****************************/
    	    	 sumVel = 0.0;
    	    	 for(int avg_i=0; avg_i<dmp->nx; avg_i++)
    	    		 sumVel += dmp->vModel[j][avg_i];

    	    	 avgVel = sumVel/dmp->nx;
				 w_v = fwp->wT[i]/avgVel;
				 w_v2 = w_v * w_v;
				 for(int l=0; l<dmp->nk; l++)
				 {
					 Kx2 = fwp->wK[l]*fwp->wK[l];
					 shift = (w_v>=fabs(fwp->wK[l]))?(cexp((sqrt(w_v2-Kx2)) * I * dmp->dz)) : (0.0+0.0*I);
					 in2[l] = out1[l] * shift;
				 } // End of Kx axis

				 fftw_execute(p2);
				 for(int m=0; m<dmp->nx; m++)
				 {
					 in1[m] = out2[m]/dmp->nk;
					 dmp->image[j][m] += crealf(in1[m]);
				 }
			 }
			 /*********End of Phase shift for reference velocity ********************/

    	     if(dmp->migType == PSPI)
    	     {
    	    	 /********* Phase shift due to reference velocity **************************/
				 for(int k=0; k<rvp->maxRefVel; k++)
				 {
					 w_v = fwp->wT[i]/rvp->refVel[j][k];
					 w_v2 = w_v * w_v;

					 for(int l=0; l<dmp->nk; l++)
					 {
						 Kx2 = fwp->wK[l]*fwp->wK[l];
						 shift = (w_v>=fabs(fwp->wK[l]))?(cexp((sqrt(w_v2-Kx2)) * I * dmp->dz)) : (0.0+0.0*I);
						 in2[l] = out1[l] * shift;
					 } // End of Kx axis

					 fftw_execute(p2);
					 for(int m=0; m<dmp->nx; m++)
						 fwp->fzData[k][m] = out2[m]/dmp->nk;
				  }
              }
    	     /*********End of Phase shift for reference velocity ***********************/

    	     /** Start performing linear interpolation ***********************************/
			 if(rvp->maxRefVel > 1 && dmp->migType == PSPI )
			 {
				 dv = (rvp->refVel[j][rvp->maxRefVel-1] - rvp->refVel[j][0])/(rvp->maxRefVel-1);
				 /*  Different velocities in a single velocity slice (dz) ***************/
				 if(dv > 0.0)
				 {
					 for(int k=0; k<dmp->nx; k++)
					 {
						 velIndx1 = (int)((dmp->vModel[j][k] - rvp->refVel[j][0])/(dv+eps));
						 velIndx2 = velIndx1+1;

						 /****** y = ( y1 (x - x0) + y0 (x1-x) ) / ( x1 - x0 ) **********/
						 if(fabs(dmp->vModel[j][k] - rvp->refVel[j][velIndx1]) > THRESHVEL && fabs(rvp->refVel[j][velIndx2] - dmp->vModel[j][k]) > THRESHVEL)
						 {
							 in1[k] = ((rvp->refVel[j][velIndx2] - dmp->vModel[j][k])*(fwp->fzData[velIndx1][k]) + (fwp->fzData[velIndx2][k]) * (dmp->vModel[j][k]-rvp->refVel[j][velIndx1]))/dv;
							 dmp->image[j][k] += crealf(in1[k]);
						 }
						 else
						 {
							 in1[k] = ((dmp->vModel[j][k] - rvp->refVel[j][velIndx1]) <= THRESHVEL)?(fwp->fzData[velIndx1][k]):(fwp->fzData[velIndx2][k]);
							 dmp->image[j][k] += crealf(in1[k]);
						 }
					 }
				 }
				 /*  Same velocity in a single velocity slice ***************************/
				 else
				 {
					 for(int k=0; k<dmp->nx; k++)
					 {
						 in1[k] = (fwp->fzData[1][k]);
						 dmp->image[j][k] += crealf(fwp->fzData[1][k]);
					 }
				 }
			  } /***** End of performing linear interpolation **********/
           } // End of z -axis
       } // End of w -axis

       void fftw_cleanup(void);
}

/* Write binary file as an Image ************************************************************/
inline void SeismicImaging::writeBinFile(std::string fileName,  int n1, int n2, float **data)
{
	int fsize;
	FILE *fp;
	fp = fopen(fileName.c_str(), "rb");
	float *temp;

	fsize = n1 * n2;
	temp = (float *) malloc(sizeof(float) * fsize);

	for(int i=0; i<n2; i++)
		for(int j=0; j<n1; j++)
			temp[i*n1+j] = data[j][i];

	fp = fopen(fileName.c_str(), "w");
	fwrite(temp,sizeof(float),fsize,fp);
	fclose(fp);
}
