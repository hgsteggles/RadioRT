/*
 * Fluid.cpp
 *
 *  Created on: 18 Dec 2014
 *      Author: harry
 */

#include "Fluid.hpp"

#include "Constants.hpp"
#include "Logger.hpp"

#include <fstream>

Fluid::Fluid(std::string filename, int ndims, const std::array<int, 3>& ncells, double sideLength)
: cells(ncells[0], ncells[1], ncells[2])
{
	std::fstream myfile(filename, std::ios_base::in);
	if (!myfile) {
		throw std::runtime_error("Fluid::Fluid: Could not input data file: " + filename + ".");
	}

	double ignore;
	std::array<int, 3> ncells_check;
	myfile >> ignore >> ncells_check[0] >> ncells_check[1] >> ncells_check[2];
	if (ncells[0] != ncells_check[0] || ncells[1] != ncells_check[1] || ncells[2] != ncells_check[2]) {
		throw std::runtime_error("Fluid::Fluid: Torch params and data have mismatched ncells.");
	}

	std::array<double, 3> x0 = std::array<double, 3>{0, 0, 0};

	dx = sideLength / ncells[0];

	for (int row = 0; row < ncells[0] * ncells[1] * ncells[2]; ++row) {
		for (int idim = 0; idim < ndims; ++idim)
			myfile >> x0[idim];
		int i = (int)(x0[0] / dx);
		int j = (int)(x0[1] / dx);
		int k = (int)(x0[2] / dx);

		if (i < 0 || j < 0 || k < 0 || i >= ncells[0] || j >= ncells[1] || k >= ncells[2]) {
			throw std::runtime_error("Fluid::Fluid: Torch params and data have mismatched deltax.");
		}

		Cell& c = cells(i, j, k);

		myfile >> c.density >> c.pressure >> c.hii;
		for (int idim = 0; idim < ndims; ++idim)
			myfile >> c.velocity[idim];
	}
	myfile.close();

	for (int i = 0; i < ncells[0]; ++i) {
		for (int j = 0; j < ncells[1]; ++j) {
			for (int k = 0; k < ncells[2]; ++k) {
				Cell& cell = cells(i, j, k);

				cell.dx = dx;
				cell.massFractions[SPECIES::H] = 1.0;
				cell.temperature = calcTemperature(cell.hii, cell.pressure,
												   cell.density,
												   cell.massFractions[SPECIES::H]);
			}
		}
	}

	updateFreeFreeCoeffs();
}

void Fluid::updateFreeFreeCoeffs() {
	double pi = Constants::PI();
	double me = Constants::electronMass();
	double e = Constants::electronCharge();
	double kb = Constants::boltzmann();
	double c = Constants::lightSpeed();
	double mp = Constants::hydrogenMass();

	double emission_ff_coeff_part = (std::pow(2.0*e*e/c, 3)/(3.0*me))
									* std::sqrt(2.0*pi/(3.0*kb*me));
	double absorption_ff_coeff_part = (4.0*std::pow(e, 6)/(3.0*me*kb*c))
									  * std::sqrt(2.0*pi/(3.0*kb*me));

	for (int i = 0; i < cells.sizex(); ++i) {
		for (int j = 0; j < cells.sizey(); ++j) {
			for (int k = 0; k < cells.sizez(); ++k) {
				Cell& cell = cells(i, j, k);

				double nhii = cell.hii * cell.massFractions[SPECIES::H] * (cell.density / mp);
				double nhe = cell.massFractions[SPECIES::HE] * (cell.density / mp) * 0.25;
				double ncno = cell.massFractions[SPECIES::CNO] * (cell.density / mp)
							  * (0.0702247191); // Number is 14.24**(-1)

				std::array<double, 4> npl = getElectrons(nhii, nhe, ncno, cell.temperature);
				double ion[] = {1.0, 2.0, 3.0, 6.0};  //assumed ionization levels
				double ne = ion[0] * npl[0] + ion[1] * npl[1] + ion[2] * npl[2] + ion[3] * npl[3];

				cell.absorption_ff_coeff = absorption_ff_coeff_part * ne / std::pow(cell.temperature, 1.5);
				cell.emission_ff_coeff = emission_ff_coeff_part * ne * std::pow(cell.temperature, -0.5);
			}
		}
	}
}

void Fluid::updateLineData(int nlevel, double freq0, double turb_broad) {
	double pi = Constants::PI();
	double me = Constants::electronMass();
	double e = Constants::electronCharge();
	double kb = Constants::boltzmann();
	double c = Constants::lightSpeed();
	double mp = Constants::hydrogenMass();
	double h = Constants::h();
	double ra = Constants::rydbergConstant();
	double eh = e / h;
	double ehc = eh / c;
	nlevel = std::min(nlevel, 1);

	double An = 64.0 * std::pow(pi, 6) * me * std::pow(ehc, 3) * std::pow(eh, 3)
				* std::pow(e, 4) / (3.0 * std::pow(nlevel, 5));
	double absorption_l_coeff = An * c * c / (8.0 * pi * freq0 * freq0);
	double Nn_coeff = nlevel * nlevel * std::pow((h * h / (2.0 * pi * me)), 1.5);
	double chi = h * c * ra / (nlevel * nlevel);

	for (int i = 0; i < cells.sizex(); ++i) {
		for (int j = 0; j < cells.sizey(); ++j) {
			for (int k = 0; k < cells.sizez(); ++k) {
				Cell& cell = cells(i, j, k);

				double nhii = cell.hii * cell.massFractions[SPECIES::H] * (cell.density / mp);
				double nhe = cell.massFractions[SPECIES::HE] * (cell.density / mp) * 0.25;
				double ncno = cell.massFractions[SPECIES::CNO] * (cell.density / mp) * (0.0702247191); // Number is 14.24**(-1)

				double tem = cell.temperature;
				std::array<double, 4> npl = getElectrons(nhii, nhe, ncno, tem);
				double ion[] = {1.0, 2.0, 3.0, 6.0};  //assumed ionization levels
				double ne = ion[0] * npl[0] + ion[1] * npl[1] + ion[2] * npl[2] + ion[3] * npl[3];

				double Nn = Nn_coeff * std::pow(kb * tem, -1.5) * nhii * ne * std::exp(-chi / (kb * tem));

				cell.absorption_l_coeff = absorption_l_coeff * Nn * (1.0 - std::exp(-h * freq0 / (kb * tem)));

				cell.sigma = std::sqrt(0.5 * (freq0 / c) * (freq0 / c)
									   * ((2 * kb * cell.temperature / mp) + turb_broad * turb_broad));
			}
		}
	}

	recl.setNLevel(nlevel);
	recl.setLineFrequency(freq0);
}

void Fluid::updateDepartureCoeffs(const std::string& filename) {

	if (filename.compare("") != 0) {
		std::ifstream file(filename, std::ios_base::in);
		if (!file)
			throw std::runtime_error("DataReader::readDataParameters: invalid input file "
									 + filename + ".");
		int nx, ny;
		file >> nx >> ny;


		std::vector<std::vector<double>> bsubn(nx, std::vector<double>(ny, 0));
		std::vector<std::vector<double>> bsubnp1(nx, std::vector<double>(ny, 0));
		std::vector<double> bsubtem(nx, 0);
		std::vector<double> bsubden(ny, 0);


		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j)
				file >> bsubtem[i] >> bsubden[j] >> bsubn[i][j] >> bsubnp1[i][j];
		}

		LogSplineData2D bsubnSpline2D = LogSplineData2D(bsubtem, bsubden, bsubn);
		LogSplineData2D bsubnp1Spline2D = LogSplineData2D(bsubtem, bsubden, bsubnp1);

		double mp = Constants::hydrogenMass();

		for (int i = 0; i < cells.sizex(); ++i) {
			for (int j = 0; j < cells.sizey(); ++j) {
				for (int k = 0; k < cells.sizez(); ++k) {
					Cell &cell = cells(i, j, k);

					double nh = cell.massFractions[SPECIES::H] * cell.density / mp;
					try {
						cell.bn = bsubnSpline2D.interpolate(cell.temperature, nh);
						cell.bnp1 = bsubnp1Spline2D.interpolate(cell.temperature, nh);
					}
					catch (std::exception& e) {
						cell.bn = 1.0;
						cell.bnp1 = 1.0;
					}
				}
			}
		}

	}
	else {
		throw std::runtime_error("departure coefficients filename cannot be empty.");
	}
}

double Fluid::getDeltaX() {
	return dx;
}

int Fluid::getNCells(int i) {
	return cells.size(i);
}

Cell const & Fluid::getCell(int i, int j, int k) {
	return cells(i, j, k);
}

Bremsstrahlung& Fluid::getBremsstrahlung() {
	return brem;
}

RecombinationLine& Fluid::getRecombinationLine() {
	return recl;
}

std::array<double, 4> Fluid::getElectrons(double nhii, double nhe, double ncno,
										  double tem) {
	std::array<double, 4> npl;
	if (tem <= 20000.0) {
		npl[0] = nhii + nhe;              // Standard CWB for Paper II
		npl[1] = ncno;
		npl[2] = 0.0;
		npl[3] = 0.0;
	}
	else if (tem <= 30000.0) {
		npl[0] = nhii + nhe;
		npl[1] = ncno;
		npl[2] = 0.0;
		npl[3] = 0.0;
	}
	else if (tem <= 1.0e6) {
		npl[0] = nhii;
		npl[1] = nhe;
		npl[2] = ncno;
		npl[3] = 0.0;
	}
	else {  // T > 1.06eK - assume CNO on average 6+ ionized
		npl[0] = nhii;
		npl[1] = nhe;
		npl[2] = 0.0;
		npl[3] = ncno;
	}
	return npl;
}

double Fluid::calcTemperature(double hii, double pre, double den,
							  double massFractionH) {
	double mu_inv = massFractionH*(hii + 1.0) + (1.0 - massFractionH)*0.25;
	return (pre/den)*(1.0/mu_inv)/Constants::specificGas();
}

void Fluid::flip() {
	cells.flip();
}


