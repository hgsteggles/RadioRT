/*
 * RecombinationLine.hpp
 *
 *  Created on: 13 Jul 2015
 *      Author: harry
 */

#ifndef RECOMBINATIONLINE_HPP_
#define RECOMBINATIONLINE_HPP_

#include "Ray.hpp"
#include "DataCube.hpp"
#include "SplineData.hpp"

class Cell;

class RecombinationLine : public Ray {
public:
	virtual void update(const Cell& cell, double ds,double freq);

	void setNLevel(int nl);
	void setLineFrequency(double freq);

private:
	int nlevel = 0;
	double freq0 = 0;
};



#endif /* RECOMBINATIONLINE_HPP_ */
