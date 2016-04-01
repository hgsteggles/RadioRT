/*
 * RecombinationLine.cpp
 *
 *  Created on: 13 Jul 2015
 *      Author: harry
 */

#include "RecombinationLine.hpp"
#include "Fluid.hpp"
#include "Constants.hpp"
#include "GaussianIntegrator.hpp"
#include "Gaunt.hpp"
#include "faddeeva/Faddeeva.hpp"

#include <cmath>
#include <iostream>
#include <fstream>

void RecombinationLine::update(const Cell& cell, double ds, double freq) {
	double c = Constants::lightSpeed();
	double pi = Constants::PI();
	double mp = Constants::hydrogenMass();
	double h = Constants::h();
	double kb = Constants::boltzmann();
	// first calculate number densities of +, ++, +++, 6+, e- (npl1, npl2, npl3, npl6, ne respectively)

	double nhii = cell.hii * cell.massFractions[SPECIES::H] * (cell.density / mp);
	double nhe = cell.massFractions[SPECIES::HE] * (cell.density / mp) * 0.25;
	double ncno = cell.massFractions[SPECIES::CNO] * (cell.density / mp) * (0.0702247191); // Number is 14.24**(-1)

	double tem = cell.temperature;

	double kbT = kb*tem;
	double exphvkt = std::exp(h*freq/kbT);
	double expmhvkt = 1.0/exphvkt;

	std::array<double, 4> npl = Fluid::getElectrons(nhii, nhe, ncno, tem);
	std::array<double, 4> ion = std::array<double, 4>{1.0, 2.0, 3.0, 6.0};  //assumed ionization levels
	std::array<double, 4> g;
	for (int n = 0; n < 4; n++)
		g[n] = gaunt(tem, ion[n], freq);

	// f-f emission and absorption coeffs (Rybicki & Lightman eqn. 5.14b + 5.19b)
	double ni_gff = (ion[0]*ion[0]*g[0]*npl[0] + ion[1]*ion[1]*g[1]*npl[1] + ion[2]*ion[2]*g[2]*npl[2] + ion[3]*ion[3]*g[3]*npl[3]);
	double emission_ff = cell.emission_ff_coeff*expmhvkt*ni_gff;
	double absorption_ff = cell.absorption_ff_coeff*ni_gff/(freq*freq);

	double sigma = cell.sigma;

	double prof = (1.0/(Constants::sqrt2pi*sigma))*std::exp(-(freq - freq0)*(freq - freq0)/(2.0*sigma*sigma));

	/*
	double faddeeva_delta = 4.7*std::pow(nlevel/100.0, 4.4)*std::pow(tem/1.0e4, 0.1)*ne;
	std::complex<double> faddeeva_z(freq - freq0, faddeeva_delta);
	faddeeva_z /= sigma*std::sqrt(2);
	double faddeeva_prof = Faddeeva::w(faddeeva_z, 1.0e-8).real()/(sigma*sqrt2pi);
	 */

	// line emission and absorption coefficients
	double absorption_l = cell.absorption_l_coeff*prof; //!< Absorption coefficient [cm-1].
	double Bv = 2.0*h*freq*freq*freq/(c*c); //!< Spectral Radiance [erg.s-1.cm-2.Hz-1.sr-1].
	Bv *= h*freq/kbT < 1.0e-5 ? kbT/(h*freq) : 1.0/(exphvkt - 1.0);
	double emission_l = Bv*absorption_l;

	// departure coefficients
	double beta = (1.0 - (cell.bnp1/cell.bn)*expmhvkt)/(1.0 - expmhvkt);
	emission_l *= cell.bnp1;
	absorption_l *= cell.bn*beta;

	double S = (emission_ff + emission_l)/(absorption_ff + absorption_l);

	double dtau = (absorption_l + absorption_ff)*ds*cell.dx;
	tau += dtau;
	intensity += (dtau > 1.0e-4) ? (1.0 - std::exp(-dtau))*(S - intensity) : dtau*(S - intensity);

	double ne = ion[0]*npl[0] + ion[1]*npl[1] + ion[2]*npl[2] + ion[3]*npl[3];
	emissionMeasure += ne*ne*ds*cell.dx;
}

void RecombinationLine::setNLevel(int nl) {
	nlevel = nl;
}

void RecombinationLine::setLineFrequency(double freq) {
	freq0 = freq;
}