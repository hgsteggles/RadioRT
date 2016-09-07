/*
 * Grid.hpp
 *
 *  Created on: 12 Dec 2014
 *      Author: harry
 */

#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>
#include <string>
#include <array>

#include "DataCube.hpp"
#include "CellGrid.hpp"
#include "Bremsstrahlung.hpp"
#include "RecombinationLine.hpp"


class Fluid {
public:
	Fluid(std::string filename, int ndims, const std::array<int, 3>& ncells, double sideLength);

	void updateFreeFreeCoeffs();
	void updateLineData(int nlevel, double freq0, double turb_broad);
	void updateDepartureCoeffs(const std::string& filename);

	double getDeltaX();
	int getNCells(int i);
	Cell const & getCell(int i, int j, int k);
	Bremsstrahlung& getBremsstrahlung();
	RecombinationLine& getRecombinationLine();

	static std::array<double, 4> getElectrons(double nhii, double nhe, double ncno, double tem);
	static double calcTemperature(double hii, double pre, double den, double massFractionH);

	void flip();

private:
	double dx = 0;
	CellGrid cells;
	Bremsstrahlung brem;
	RecombinationLine recl;
};



#endif /* GRID_HPP_ */
