/*
 * RayTracer.hpp
 *
 *  Created on: 12 Dec 2014
 *      Author: harry
 */


#ifndef RAYTRACER_HPP_
#define RAYTRACER_HPP_

#include "Ray.hpp"
#include "DataCube.hpp"
#include "Bremsstrahlung.hpp"
#include "RecombinationLine.hpp"

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <array>
#include <memory>
#include <vector>

class Fluid;

/*
class RayTracer {
public:
	std::vector<std::vector<double>> image;
	std::vector<std::vector<double>> tau;
	double inten_ff;
	double fac;

	RayTracer(int k, double theta, double phi, int xpixels, int ypixels, double sampling, double freq, const Fluid& fluid);
	RayTracer(double theta, int xpixels, int ypixels, double sampling, double freq, const Fluid& fluid);
};
*/

class RayTracerData {
public:
	DataCube fluxFF;
	DataCube tauFF;
	DataCube emissionMeasureFF;
	double intensityFF = 0;

	DataCube fluxRL;
	DataCube tauRL;
	DataCube emissionMeasureRL;
	double intensityRL = 0;

	double fac = 0;
};

class TrigData {
public:
	TrigData(double theta);
	double sin() const;
	double cos() const;
	double tan() const;
	double ang() const;
private:
	double cos_th, sin_th, tan_th, angle;
};

class RayTracer {
public:
	double calcPixelSizeModification();
	static std::array<int, 2> calcPixelNumberAxi(int rcells, int zcells, double dx, double desiredPixelSize);
	double calcPixelSizeModificationAxi(int rcells, int zcells);

	void setImageDimensions(int xpixels, int ypixels);
	void setSampling(int s);

	void setFrequencies(const std::vector<double>& freqs);
	void setIntegratingRL(bool on);
	void setIntegratingFF(bool on);
	void setDopplerShifted(bool dopplerShifted);
	void setLineOfSightVelocity(double vLOS);
	void setDopplerIntMinPhiInc(double doppShiftPhiIncr);

	RayTracerData rayTrace3D(Fluid& fluid, int k, double theta, double phi);
	RayTracerData rayTraceAxiSymm(Fluid& fluid, double theta);

private:
	std::array<int, 2> pixels;
	int sampling = 4;

	std::vector<double> frequencies;
	bool integratingRL = false;
	bool integratingFF = false;
	bool isDopplerShifted = false;
	double vLOS = 0;
	double doppShiftPhiIncr = 0.01; //!< Increment in azimuthal angle when calculating doppler shifts across cell.
};

class DopplerShift {
public:
	DopplerShift(double frequency, double vLOS, double ur, double uz, double x, const TrigData& trig_th, double phi1, double phi2, double path_length, double min_angle);

	std::vector<double> shifts;
	std::vector<double> deltas;
};


#endif /* RAYTRACER_HPP_ */
