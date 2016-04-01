/*
 * Bremsstrahlung.hpp
 *
 *  Created on: 13 Jul 2015
 *      Author: harry
 */

#ifndef BREMSSTRAHLUNG_HPP_
#define BREMSSTRAHLUNG_HPP_

#include "Ray.hpp"
#include "DataCube.hpp"

class Cell;

class Bremsstrahlung : public Ray {
public:
	virtual void update(const Cell&, double ds, double freq);

private:
};



#endif /* BREMSSTRAHLUNG_HPP_ */
