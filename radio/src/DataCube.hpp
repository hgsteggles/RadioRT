/*
 * DataCube.hpp
 *
 *  Created on: 15 Jul 2015
 *      Author: harry
 */

#ifndef DATACUBE_HPP_
#define DATACUBE_HPP_

#include <vector>
#include <cmath>
#include <iostream>

class DataCube {
public:
	DataCube()
	: DataCube(1, 1, 1)
	{

	}

	DataCube(const unsigned int nx, const unsigned int ny, const unsigned int nz)
	: nx(nx)
	, ny(ny)
	, nz(nz)
	, vector3D(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)))
	{

	}

	unsigned int sizex() const {
		return nx;
	}
	unsigned int sizey() const {
		return ny;
	}
	unsigned int sizez() const {
		return nz;
	}

	double& operator()(const unsigned int ix, const unsigned int iy, const unsigned int iz) {
		return vector3D[ix][iy][iz];
	}

	const double& operator()(const unsigned int ix, const unsigned int iy, const unsigned int iz) const {
		return vector3D[ix][iy][iz];
	}

	DataCube operator* (double x) {
		DataCube cube = *this;

		for (unsigned int i = 0; i < cube.sizex(); ++i) {
			for (unsigned int j = 0; j < cube.sizey(); ++j) {
				for (unsigned int k = 0; k < cube.sizez(); ++k) {
					cube(i, j, k) *= x;
				}
			}
		}

		return cube;
	}

	DataCube& mul(double x) {
		for (unsigned int i = 0; i < sizex(); ++i) {
			for (unsigned int j = 0; j < sizey(); ++j) {
				for (unsigned int k = 0; k < sizez(); ++k) {
					vector3D[i][j][k] *= x;
				}
			}
		}

		return *this;
	}

private:
	unsigned int nx = 0, ny = 0, nz = 0;
	std::vector<std::vector<std::vector<double>>> vector3D;
};

#endif /* DATACUBE_HPP_ */
