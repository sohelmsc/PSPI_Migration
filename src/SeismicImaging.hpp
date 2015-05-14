/*  =================================================================================
  *  Name       : SeismicImaging.hpp
 *  Author      : Mafijul Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Base code of this software package. Contains all the basic methods entailed to generate image.
 *  			: This is a generalised migration code which can include different migration methods. It can accept
 *  			: any number of new method declaration.
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef SEISMICIMAGING_HPP_
#define SEISMICIMAGING_HPP_

#include "DataModelParams.hpp"
#include "FreqWavenumberParams.hpp"
#include "RefVelParams.hpp"

class SeismicImaging{

	/**
	 * Public methods:
	 * Indeed, I encapsulated different types of migration in different public methods.
	 * The naming of those methods are self-explanatory.
	 **/
	public:
			void Imaging(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp);
			//void stoltImaging(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp);
			//void splitStepImaging(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp);

	/**
	 * Access modifier of these methods are made private to prohibit the
	 * execution of these functions from anywhere except the object of this class
	 **/
	private:
        /* Read binary file. n1--> Fast direction (time), n2-> slow direction (space), data->data stored in this variable */
		int readBinFile(std::string fileName,  int n1, int n2, float **data);
		float ** alloc2DPointer(int n1, int n2, float **data);
		float complex ** alloc2DComPointer(int n1, int n2, float complex **data);
		void computeRefVel(RefVelParams* rvp, DataModelParams* dmp);
		void tXtoFx(DataModelParams* dmp, FreqWavenumberParams* fwp);
		void computeFreqWavenumber(FreqWavenumberParams* fwp, DataModelParams* dmp);
		void downwardExtrapolation(DataModelParams* dmp, FreqWavenumberParams* fwp, RefVelParams* rvp1);
		void fXtoFk(DataModelParams* dmp, FreqWavenumberParams* fwp);
		void writeBinFile(std::string fileName,  int n1, int n2, float **data);
};
#endif /* SEISMICIMAGING_HPP_ */
