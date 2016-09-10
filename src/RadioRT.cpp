/*
 * RadioRT.cpp
 *
 *  Created on: 8 Jan 2015
 *      Author: harry
 */

#include "RadioRT.hpp"
#include "Timer.hpp"
#include "Constants.hpp"
#include "Converter.hpp"
#include "Fluid.hpp"
#include "Bremsstrahlung.hpp"
#include "RecombinationLine.hpp"
#include "RayTracer.hpp"
#include "WriteFITS2D.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <exception>
#include <array>
#include <memory>

RadioRT::RadioRT(const RadioRT_Parameters& p)
: params(p)
{

}

void RadioRT::run() {
	Timer timer;
	timer.start();

	Fluid fluid(params.inputfile, params.ndims, params.ncells, params.sideLength);
	if (params.integratingRL) {
		fluid.updateLineData(params.nlevel, params.frequency, params.turb_broad);
		fluid.updateDepartureCoeffs("config/bsubn.txt");
	}

	if (params.theta < -90.0 || params.theta > 90.0)
		throw std::runtime_error("RadioRT::run: theta is out of range (0 <= th <= 180).");
	if (params.phi < 0.0 || params.phi > 360.0)
		throw std::runtime_error("RadioRT::run: phi is out of range (0 <= ph <= 360).");

	if (params.theta < 90.0)
		fluid.flip();

	// Set the number of pixels in the ray-traced image (npixels[0] and npixels[1])
	std::array<int, 2> npixels;
	npixels[0] = params.ncells[0] * params.resolutionScale;
	npixels[1] = params.ncells[1] * params.resolutionScale;
	if (params.geometry == "cylindrical") {
		npixels[0] *= 2;
	}

	//convert theta and phi to radians
	double theta_rad = Converter::DEG_2_RAD(params.theta);
	double phi_rad = Converter::DEG_2_RAD(params.phi);

	// Calculate frequencies.
	std::vector<double> freqs;
	double channelWidth = params.bandwidth / params.nchannels;
	for (int i = 0; i < params.nchannels; ++i)
		freqs.push_back(params.frequency + (i + 0.5) * channelWidth - (params.bandwidth / 2.0));

	// make the image on the plane of the sky and return it in the array image(npixels[0],npixels[1]) (units are erg/s/cm^2/Hz/ster).
	// The opacity is returned in the array tau(npixels[0],npixels[1])
	RayTracer raytracer;
	raytracer.setIntegratingFF(params.integratingFF);
	raytracer.setIntegratingRL(params.integratingRL);
	raytracer.setImageDimensions(npixels[0], npixels[1]);
	raytracer.setSampling(params.sampling);
	raytracer.setFrequencies(freqs);
	raytracer.setDopplerShifted(params.dopplerShifted);
	raytracer.setLineOfSightVelocity(params.vLOS);
	raytracer.setDopplerIntMinPhiInc(params.doppShiftPhiIncr);

	RayTracerData data;
	if (params.geometry == "cylindrical")
		data = raytracer.rayTraceAxiSymm(fluid, std::abs(theta_rad));
	else
		data = raytracer.rayTrace3D(fluid, std::abs(theta_rad), phi_rad);

	if (params.theta < 90.0)
		data.flip();

	// Scale to correct units and work out total flux
	double pixsize = data.fac * fluid.getDeltaX(); // pixel size [cm].
	double pixrad = pixsize / (Units::milli() * params.dist * Converter::PC_2_CM(1)); // pixel size [radians].
	double pixdeg = Converter::RAD_2_DEG(pixrad); // pixel size [degrees].
	double pixarcm = pixdeg * 60.0; // pixel size [arcminutes].
	double pixarcs = pixarcm * 60.0; // pixel size [arcseconds].
	double pixmas = Units::milli() * pixarcs; // pixel size [milliarcseconds].

	data.intensityFF *= Units::milli() * Converter::CGS_2_JY(1) * pixrad * pixrad; //mJy
	data.intensityRL *= Units::milli() * Converter::CGS_2_JY(1) * pixrad * pixrad;
	data.fluxFF.mul(Units::milli() * Converter::CGS_2_JY(1) * pixrad * pixrad);
	data.fluxRL.mul(Units::milli() * Converter::CGS_2_JY(1) * pixrad * pixrad);

	data.emissionMeasureFF.mul(Converter::CM_2_PC(1.0));
	data.emissionMeasureRL.mul(Converter::CM_2_PC(1.0));

	if (params.integratingFF) {
		WriteFITS::wfits(params.outputDirectory + "odepth_ff", "", data.tauFF, params, pixdeg, pixsize, fluid.getDeltaX());
		WriteFITS::wfits(params.outputDirectory + "emeasure_ff", "cm**6/pc", data.emissionMeasureFF, params, pixdeg, pixsize, fluid.getDeltaX());
		WriteFITS::wfits(params.outputDirectory + "intensity_pixel_ff", "mJy/pixel", data.fluxFF, params, pixdeg, pixsize, fluid.getDeltaX(), data.intensityFF);
	}
	if (params.integratingRL) {
		WriteFITS::wfits(params.outputDirectory + "odepth_rl", "", data.tauRL, params, pixdeg, pixsize, fluid.getDeltaX());
		WriteFITS::wfits(params.outputDirectory + "emeasure_rl", "cm**6/pc", data.emissionMeasureRL, params, pixdeg, pixsize, fluid.getDeltaX());
		WriteFITS::wfits(params.outputDirectory + "intensity_pixel_rl", "mJy/pixel", data.fluxRL, params, pixdeg, pixsize, fluid.getDeltaX(), data.intensityRL);
	}
}

