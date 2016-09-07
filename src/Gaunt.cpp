/*
 * Gaunt.cpp
 *
 *  Created on: 18 Dec 2014
 *      Author: harry
 */

#include "Gaunt.hpp"

#include <cmath>
#include <stdexcept>
#include <iostream>

double gaunt(double tem, double charge, double freq) {
	// Chebyshev coefficients
	double d[8][11] = {
		{8.986940175e00, -8.006936989e-1, -3.781305103e-1,
		1.877213231e-2,  7.300158392e-2, -1.744671550e-3,
		-1.707268366e-2,  2.567331664e-4,  4.098322531e-3,
		3.837562402e-5, -8.491991820e-4},
		{-4.009515855e00,  9.466021705e-1,  1.102726322e-1,
		-1.004885705e-1,  3.576785497e-3,  2.864013856e-2,
		-4.694254776e-3, -9.155339970e-3,  1.635218463e-3,
		2.938325230e-3, -3.615327726e-4},
		{8.808871266e-1,  9.043402532e-2, -1.543619180e-2,
		-5.483366378e-2, -4.545307025e-3,  1.903394837e-2,
		1.311691517e-3, -6.997479192e-3, -5.918883504e-4,
		2.393747064e-3,  3.148015257e-4},
		{2.640245111e-2, -9.608451450e-2,  8.310561114e-3,
		-4.520154409e-3, -1.017965604e-2,  7.091074494e-3,
		5.316703136e-3, -3.571518641e-3, -2.333091048e-3,
		1.328839809e-3,  8.909207650e-4},
		{-4.580645915e-2, -1.885629865e-2,  2.179620525e-2,
		8.366530426e-3, -9.530211924e-3, -9.668371391e-4,
		5.178193095e-3, -2.096101038e-4, -2.484138313e-3,
		9.135013312e-5,  9.869737522e-4},
		{-3.568055702e-3,  1.050313890e-2,  4.259726289e-3,
		3.700273930e-3, -3.450186162e-3, -2.999107465e-3,
		2.451228935e-3,  1.553822487e-3, -1.359996060e-3,
		-7.137252303e-4,  6.134671184e-4},
		{2.827798067e-3,  2.800889961e-3, -4.181588794e-3,
		6.889320423e-4,  1.040482914e-3, -1.820642230e-3,
		-2.277321615e-5,  1.509584686e-3, -5.371426147e-5,
		-7.656848158e-4,  1.068883394e-4},
		{3.365860195e-4, -1.078209202e-3, -1.770208330e-3,
		9.460313195e-5,  1.407073544e-3, -3.874082085e-4,
		-8.182359057e-4,  6.212627837e-4,  5.553549563e-4,
		-3.504683798e-4, -2.046080100e-4}
	};

	// calculate gamma2 gamma and u from z, nu and T

	double gamma2 = 157833.0*charge*charge/tem;   // introduced gamma2, used a few times
	double gamma = std::sqrt(gamma2);
	double u = 4.797978e-11*freq/tem;

	// check if parameters are in range for Hummer's fit
	if (gamma2 >= 1.0e-3 && gamma2 <= 1.0e3 && u >= 1.0e-4 && u <= 31.622777) {
		// normalise gamma and u onto -1, 1 and evaluate the
		// Chebyshev polynomials
		double xu = (2.0*std::log10(u) + 2.5)*0.18181818181818181818;
		double xg = std::log10(gamma2)*0.333333333333333333333;


		double cj[8];
		for (int j = 0; j < 8; ++j) {
			double dcol[11];
			for (int i = 0; i < 11; ++i)
				dcol[i] = double(d[j][i]);
			cj[j] = chebev(xg, dcol, 11);
		}
		return chebev(xu,cj,8);
	}


	// not in range for Hummer's fit - check Scheuer's approximation is ok.
	if (u < 1.0e-4 && gamma >= 1.0) {
		// use Scheuer's (1960) long-wavelength approximation (see Hummer).
		return -0.55133*(std::log(gamma) + std::log(u) + 0.056745);
	}

	// not in range for Scheuer's fit, try Elwert's high-energy approx (see Hummer).
	if (u < 1.0e-4 && gamma < 1.0) {
		// use Elwert's (1954) approximation (see Hummer).
		return 0.55133*(0.80888 - std::log(u));
	}

	// give up and provide a warning!

	// The Gaunt factor gets set to 1 here because in this
	// particular problem we want a low-temperature sphere
	// representing a star to be optically thick.

	//g = 1.0; //this is not too bad an approximation in many cases if
	//we are not too worried about complete accuracy - see
	//eg Hummer (1988) and Grant (1958)

	//use the asymptotic formula of Menzel & Perkeris (see Hummer 1988)
	/*
	std::cout << "Warning: gaunt factor not calculated properly." << std::endl;
	std::cout << "u = " << u << std::endl;
	std::cout << "f = " << freq << std::endl;
	std::cout << "T = " << tem << std::endl;
	std::cout << "gamma2 = " << gamma2 << std::endl;
	exit(0);
	*/
	return 1.0;
	return 1.0 + 0.1728*std::pow(u/gamma2,0.3333)*(1.0+2.0/u) - 0.0496*(u/gamma2,0.6667)*(1.0 + (0.6667/u) + 1.3333/(u*u));
}

double chebev(double y, double c[], int m){
	// Numerical Recipes Chebyshev polynomial evaluation routine
	if (y < -1.0 || y > 1.0)
		throw std::runtime_error("Gaunt::chebev: y not in range.");

	double d = 0.0;
	double dd = 0.;
	double y2 = 2.0*y;
	for (int j = m-1; j >= 1; --j){
		double sv = d;
		d = y2*d - dd + c[j];
		dd = sv;
	}
	double cheb = y*d - dd + 0.5*c[0];

	return cheb;
}
