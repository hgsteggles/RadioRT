/*
 * Ray.cpp
 *
 *  Created on: 14 Jul 2015
 *      Author: harry
 */

#include "Ray.hpp"

Ray::~Ray()
{

}

double Ray::getTau() {
	return tau;
}

double Ray::getIntensity() {
	return intensity;
}

double Ray::getEmissionMeasure() {
	return emissionMeasure;
}

void Ray::reset() {
	intensity = 0;
	tau = 0;
	emissionMeasure = 0;
}


