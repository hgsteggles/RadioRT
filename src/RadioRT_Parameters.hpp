/*
 * RadioRT_Parametes.hpp
 *
 *  Created on: 8 Jan 2015
 *      Author: harry
 */

#ifndef RADIORT_PARAMETES_HPP_
#define RADIORT_PARAMETES_HPP_

#include <string>
#include <array>

// Angle of los to model grid is specified by the angles theta and phi [degrees].
//
// theta =  angle of line of sight to the z-axis of the model grid (theta = 0 is the N-pole, theta=90 is in the equatorial plane,
//          theta = 180 is the S-pole)
// phi   =  angle of projected line of sight on the x-y plane to the x-axis of the model grid (phi = 0 or 360 is viewing along
//          the x-axis)
// NOTE: if theta = 0, and phi = 0, and the orbital plane lies in the xy plane, then the image on the sky will be as viewed from
//       above the orbital plane, with an orientation such that the hydro x-axis is vertical, with the positive hydro y-axis
//       increasing to the right

class RadioRT_Parameters {
public:
	// Torch Params.
	std::string torchParamFilename = "";
	int resolutionScale = 1;
	int ndims = 0;
	std::array<int, 3> ncells = std::array<int, 3>{0, 0, 0};
	std::string geometry = "cartesian";
	double sideLength = 0;

	// Data params.
	double dist = 0; //!< Distance to source [kpc].
	double rightAscension = 0; //!< Right ascension of source.
	double declination = 0; //!< Left ascension of source.

	// Recombination line params.
	double turb_broad = 0; // Turbulent broadening [cm2.s-2].
	bool dopplerShifted = false;
	double doppShiftPhiIncr = 10; //!< Increment in azimuthal angle when calculating doppler shifts across cell.
	//double channelProfIntegEpsilon = 1e-6; //!< Error tolerance in integration across channel profile.
	//int channelProfIntegMaxRecursion = 20; //!< Max no. of recursion levels in integration across channel profile.

	// Ray tracing params.
	double theta = 0; //!< Angle LOS makes with z axis [degrees].
	double phi = 0; //!< Angle LOS makes with x axis when projected on xy plane [degrees].
	double frequency = 0; //!< Observing frequency [Hz].
	double bandwidth = 0;
	int nchannels = 1;
	int stokes = 1;
	int nlevel = 0;
	double vLOS = 0; //!< Line of sight velocity of source away from observer [cm.s-1].

	bool integratingFF = false;
	bool integratingRL = false;
	double sampling = 0; //!< Sampling of ray-tracing. Artifacts disappear when pathlength ~0.25 of edge length (sampling = 4).
	double pixelResolutionArcseconds = -10;

	std::string inputfile = ""; //!< Input file name.
	std::string outputDirectory = "tmp"; //!< Output directory name.
};



#endif /* RADIORT_PARAMETES_HPP_ */
