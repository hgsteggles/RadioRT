/*
 * GaussianIntegrator.hpp
 *
 *  Created on: 27 Jul 2015
 *      Author: harry
 */

#ifndef GAUSSIANINTEGRATOR_HPP_
#define GAUSSIANINTEGRATOR_HPP_


/** Integrates exp(-([x - x0]^2)/(2*sigma^2))/(2*PI*sigma^2)^0.5
 *
 */
class GaussianIntegrator {
public:
	GaussianIntegrator(double x0, double sigma2, double epsilon, double maxRecursionDepth);
	double integrate(double a, double b);
	double integrateAux(double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
	double gaussian(double x);
private:
	double x0;
	double sigma2;
	double epsilon;
	double maxRecursionDepth;
	double PI = 3.14159265358979323846;
};



#endif /* GAUSSIANINTEGRATOR_HPP_ */
