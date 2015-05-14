/*  =================================================================================
 *  Name        : DataModelParams.hpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Store seismic data and model parameters
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef DATAMODELPARAMS_HPP_
#define DATAMODELPARAMS_HPP_

#include <iostream>
using std::cout;
using std::ostream;
using std::endl;


/**
 * This data structure holds the dimension and grid size of model, and seismic trace constant parameters
 */

enum migrationType{PSPI, GAZDAG, STOLT, SPLITSTEP }; /// Different migration types

class DataModelParams{
public:
	// These parameters depend on the seismic data and velocity model
	int nx;	        			/// Number of grid points in x (horizontal) direction
	int nz;	        			/// Number of grid points in z (vertical) direction
	int dx;          			/// Sampling interval in x (horizontal) direction
	int dz;          			/// Sampling interval in z (vertical) direction
	int nk;          			/// Number of sampling points in wave-number axis
	int nw;          			/// Number of sampling points in temporal frequency axis
	int nt;			 		/// Number of time samples in a seismic trace
	float dt;			 	/// Sampling rate in seismic data
	int nTraces;     			/// Number of seismic traces in 2-D data
	float **data;          			/// 2-D Seismic data
	float **vModel;        			/// Velocity model
	float **image;   	  		/// Final seismic image
	std::string dataFileName;   		/// Absolute path of the seismic data file
	std::string velFileName;   		/// Absolute path of the velocity model file
	std::string imageFileName;  		/// Absolute path of the image file to be written
	migrationType migType;        /// Assign migration type

	/// clear configuration variables
	/* Delete seismic data, image, and velocity model from the dmp object */

	void reset()
	{
		delete vModel;
		delete image;
		delete data;
	}
};

#endif /* VELOCITYMODELPARAMS_HPP_ */
