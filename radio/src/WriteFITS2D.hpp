/*
 * WriteFITS2D.hpp
 *
 *  Created on: 18 Dec 2014
 *      Author: harry
 */

#ifndef WRITEFITS2D_HPP_
#define WRITEFITS2D_HPP_

#include "DataCube.hpp"
#include "RadioRT_Parameters.hpp"

#include <string>

namespace WriteFITS {
	void printFITS_Error(int status);
	void wfits(std::string filename, std::string units, const DataCube &cube, RadioRT_Parameters params, double pixdeg,
			   double pixsize, double dx, double total = -10);

}

#endif /* WRITEFITS2D_HPP_ */
