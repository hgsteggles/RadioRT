//
// Created by harry on 17/02/16.
//

#ifndef RADIO_CELLGRID_HPP
#define RADIO_CELLGRID_HPP

#include <vector>
#include <algorithm>

#include "Cell.hpp"

class CellGrid {
public:
	CellGrid(const unsigned int nx, const unsigned int ny, const unsigned int nz)
	: nx(nx)
	, ny(ny)
	, nz(nz)
	, vector3D(nx, std::vector<std::vector<Cell>>(ny, std::vector<Cell>(nz)))
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

	unsigned int size(unsigned int i) {
		if (i == 0)
			return sizex();
		else if (i == 1)
			return sizey();
		else if (i == 2)
			return sizez();
		else
			return 0;
	}

	Cell& operator()(const unsigned int ix, const unsigned int iy, const unsigned int iz) {
		return vector3D[ix][iy][iz];
	}

	const Cell& operator()(const unsigned int ix, const unsigned int iy, const unsigned int iz) const {
		return vector3D[ix][iy][iz];
	}

	void flip() {
		for (int i = 0; i < sizex(); ++i)
			std::reverse(vector3D[i].begin(), vector3D[i].end());
	}

private:
	unsigned int nx = 0, ny = 0, nz = 0;
	std::vector<std::vector<std::vector<Cell>>> vector3D;
};


#endif //RADIO_CELLGRID_HPP
