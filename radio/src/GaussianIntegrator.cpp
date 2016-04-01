/*
 * GaussianIntegrator.cpp
 *
 *  Created on: 27 Jul 2015
 *      Author: harry
 */

#include "GaussianIntegrator.hpp"

#include <cmath>

GaussianIntegrator::GaussianIntegrator(double x0, double sigma2, double epsilon, double maxRecursionDepth)
: x0(x0)
, sigma2(sigma2)
, epsilon(epsilon)
, maxRecursionDepth(maxRecursionDepth)
{

}

double GaussianIntegrator::integrate(double a, double b) {
	double c = (a + b)/2, h = b - a;
	double fa = gaussian(a);
	double fb = gaussian(b);
	double fc = gaussian(c);
	double S = (h/6)*(fa + 4*fc + fb);
	return integrateAux(a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

double GaussianIntegrator::integrateAux(double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom) {
	double c = (a + b)/2.0, h = b - a;
	double d = (a + c)/2.0, e = (c + b)/2.0;
	double fd = gaussian(d);
	double fe = gaussian(e);
	double Sleft = (h/12.0)*(fa + 4.0*fd + fc);
	double Sright = (h/12.0)*(fc + 4.0*fe + fb);
	double S2 = Sleft + Sright;
	if (bottom <= 0 || std::fabs(S2 - S) <= 15.0*epsilon)
		return S2 + (S2 - S)/15.0;
	return integrateAux(a, c, epsilon/2.0, Sleft, fa, fc, fd, bottom-1) +
			integrateAux(c, b, epsilon/2.0, Sright, fc, fb, fe, bottom-1);
}

double GaussianIntegrator::gaussian(double x) {
	double result = 1.0/std::sqrt(2.0*PI*sigma2);
	result *= std::exp(-(x - x0)*(x - x0)/(2.0*sigma2));
	return result;
}


