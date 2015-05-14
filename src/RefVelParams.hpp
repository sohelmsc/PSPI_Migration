/*  =================================================================================
  *  Name        : RefVelParams.hpp
 *  Author      : Mafijul Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Store reference velocities information for PSPI
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef REFVELPARAMS_HPP_
#define REFVELPARAMS_HPP_

/**
 * This data structure holds the parameters for reference velocity information
 */

class RefVelParams{
public:
	int 	maxRefVel;     /// Maximum number of reference velocity
	int 	*nRefVel;      /// Number of reference velocities in each vertical slice
	float  **refVel;      /// Reference velocities in each vertical slice

	/******** clear configuration variables **********/
	void reset()
	{
		delete refVel;
	}
};



#endif /* REFVELPARAMS_HPP_ */
