/** Provides Constants and EnumParser classes.
 * @file Constants.hpp
 *
 *  @author "Harrison Steggles"
 *
 *  @date 24/11/2014 - the first version.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "Converter.hpp"

#include <string>
#include <map>
#include <stdexcept>
#include <cmath>

/**
 * @class Constants
 *
 * @brief Contains all constants that may be needed at any place in the code. Also has a Converter handy for unit conversions.
 *
 * @version 0.8, 24/11/2014
 */
class Constants {
public:
	Constants() = default;
	void initialise_MLT(double mass, double length, double time);
	void initialise_DPT(double density, double pressure, double time);

	Converter converter;

	//double specificGasConstant = 0; //!< Specific Gas Constant.
	//double rydbergEnergy = 0; //!< Rydberg Energy.
	//double boltzmannConst = 0; //!< Boltzmann's Constant.
	//double dustExtinctionCrossSection = 0; //!< Dust Extinction Cross Section (Baldwin et. al. 1991) [cm2 H-1].
	//double hydrogenMass = 0; //!< Mass of hydrogen [g].
	//double alphaB = 0; //!< Radiative recombination coefficient of ionised hydrogen [cm^3 s^-1].

	static double PI();
	static double lightSpeed();
	static double boltzmann();
	static double electronMass();
	static double electronCharge();
	static double rydbergEnergy();
	static double rydbergConstant();
	static double h();
	static double hbar();
	static double hydrogenMass();
	static double specificGas();

	static const double sqrt2pi;

private:
	void initialise();
};

#endif // CONSTANTS_HPP_
