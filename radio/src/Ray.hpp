/*
 * Ray.hpp
 *
 *  Created on: 13 Jul 2015
 *      Author: harry
 */

#ifndef RAY_HPP_
#define RAY_HPP_

#include <vector>
#include <string>

class Cell;

class Ray {
public:
	virtual ~Ray();
	virtual void update(const Cell&, double ds, double freq) = 0;
	void reset();
	double getTau();
	double getIntensity();
	double getEmissionMeasure();

protected:
	double tau = 0;
	double intensity = 0;
	double emissionMeasure = 0;
};



#endif /* RAY_HPP_ */
