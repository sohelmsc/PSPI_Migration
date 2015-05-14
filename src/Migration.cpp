/*  =================================================================================
 *   Name        : Migration.cpp (Phase Shift Migration Plus Interpolation and others)
 *  Author      : Mafijul Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Generate the migrated image of the velocity model using PSPI method
 *  Date        : Feb 18, 2014 - implementation for PSPI poststack method
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */


#include "DataModelParams.hpp"
#include "FreqWavenumberParams.hpp"
#include "RefVelParams.hpp"
#include "SeismicImaging.hpp"

int main()
{
   clock_t start, finish;

   /** Creating the objects of the dependence classes ************************************/
   DataModelParams dmp;
   FreqWavenumberParams fwp;
   RefVelParams rvp;

   /** Initialisation of data and model parameters. **************************************/
   dmp.dx = 10;
   dmp.dz = 5;
   dmp.nTraces = 700;
   dmp.nx = 700;
   dmp.nz = 600;
   dmp.nt = 901;
   dmp.nw = pow(2,ceil(log10(dmp.nt)/log10(2)));
   dmp.nk = pow(2,ceil(log10(dmp.nx)/log10(2)));
   dmp.dt = 0.004;
   dmp.dataFileName = "/home/entropy/workspace/PSPI_PS/src/data.bin";
   dmp.velFileName = "/home/entropy/workspace/PSPI_PS/src/vel.bin";
   dmp.imageFileName = "image_mean.bin"; //"/home/entropy/workspace/PSPI_PS/src/image.bin";
   //dmp.migType = PSPI;
   dmp.migType = GAZDAG;

  /* To embed other post-stack migration methods. Not implemented yet though.
   *
	 dmp.migType = STOLT;
	 dmp.migType = SPLITSTEP;
   *
   */

   /** Initialisation of frequency parameters. *******************************************/
   fwp.minF = 1;
   fwp.maxF = 40;

   /** Initialisation of the number of reference velocities. Should be greater than 1 for
    * PSPI method */
   rvp.maxRefVel=6;

   /** Executing migration to generate the seismic image ****************************/
   start = clock();
   SeismicImaging si;
   si.Imaging(&dmp, &fwp, &rvp);
   finish = clock();
   std::cout << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds " <<std::endl;

   /** Clearing the memory space **********************************************************/
   dmp.reset();
   fwp.reset();
   rvp.reset();
   std::cout << "Cleared the memory space." <<std::endl;

 return 0;
}
