#include "Constants.hpp"

void Constants::initialise() {
	//hydrogenMass = converter.toCodeUnits(1.674e-24, 1, 0, 0);
	//specificGasConstant = converter.toCodeUnits(8.314462e7, 0, 2, -2);
	//boltzmannConst = converter.toCodeUnits(1.3806488e-16, 1, 2, -2);
	//rydbergEnergy = converter.toCodeUnits(converter.EV_2_ERGS(13.6), 1, 2, -2);
	//dustExtinctionCrossSection = converter.toCodeUnits(5.0e-22, 0, 2, 0);
	//alphaB = converter.toCodeUnits(2.59e-13, 0, 3, -1);
}

void Constants::initialise_MLT(double mass, double length, double time) {
	converter.set_mass_length_time(mass, length, time);
	initialise();
}

void Constants::initialise_DPT(double density, double pressure, double time) {
	converter.set_rho_pressure_time(density, pressure, time);
	initialise();
}

double Constants::PI() {
	return 3.14159265358979323846;
}

double Constants::lightSpeed() {
	return 2.99792458e10;
}

double Constants::boltzmann() {
	return 1.380658e-16;
}

double Constants::electronMass() {
	return 9.1093898e-28;
}

double Constants::electronCharge() {
	return 4.8032068e-10;
}

double Constants::rydbergEnergy() {
	return 2.1798741e-11;
}

double Constants::rydbergConstant() {
	return 1.096776e5;
}

double Constants::h() {
	return 6.6260755e-27;
}

double Constants::hbar() {
	return h()/(2.0*PI());
}

double Constants::hydrogenMass() {
	return 1.674e-24;
}

double Constants::specificGas() {
	return 8.314462e7;
}

const double Constants::sqrt2pi = std::sqrt(2.0*Constants::PI());
