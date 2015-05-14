/*  =================================================================================
 *  Name        : FreqWavenumberParams.hpp
 *  Author      : Mafijul Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Store temporal and spatial frequency information
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef FREQWAVENUMBERPARAMS_HPP_
#define FREQWAVENUMBERPARAMS_HPP_

#include <complex.h>
#include <cmath>
#include <time.h>


/**
 * This data structure holds the parameters for temporal frequency and spatial Wavenumber axes
 * and Fourier transformed seismic data
 **/

class FreqWavenumberParams{
public:
	// These parameters depend on the seismic data and velocity model
	float *wT;	    			/// Temporal frequency axis
	float *wK;	    			/// Spatial wave-number axis
	float complex **fxData;	    		/// Temporal frequency data
	float complex **fkData;	    		/// Spatial wave-number data
	float complex **fzData;	    		/// Shifted f-k data
	int maxF;     				/// Maximum temporal frequency
	int minF;     				/// Minimum temporal frequency

	/******** clear configuration variables **********/
	void reset()
	{
		delete wT;
		delete wK;
		delete fxData;
		delete fzData;
	}
};

#endif /* FREQWAVENUMBERPARAMS_HPP_ */
