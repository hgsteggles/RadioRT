/** Provides the Converter class.
 * @file Converter.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

/**
 * @class Converter
 *
 * @brief Converts values to/from code units from/to cgs real units.
 *
 * @version 0.8, 24/11/2014
 */
class Converter {
public:
	Converter();
	void set_mass_length_time(const double mass, const double length, const double time);
	void set_rho_pressure_time(const double rho, const double pressure, const double time);
	double toCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const;
	double fromCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const;

	static double JY_2_CGS(double val_in_jy);
	static double CGS_2_JY(double val_in_cgs);
	static double EV_2_ERGS(double val_in_ev);
	static double ERG_2_EV(double val_in_erg);
	static double YR_2_S(double val_in_yr);
	static double S_2_YR(double val_in_s);
	static double PC_2_CM(double val_in_pc);
	static double CM_2_PC(double val_in_cm);
	static double MO_2_G(double val_in_mo);
	static double DEG_2_RAD(double val_in_deg);
	static double RAD_2_DEG(double val_in_rad);

	void printInfo();
private:
	double M, L, T;
	double convertCodeUnits(const double val, const double mass_index, const double length_index, const double time_index, const bool& from) const;
};

class Units {
public:
	static double yotta();
	static double zetta();
	static double exa();
	static double peta();
	static double tera();
	static double giga();
	static double mega();
	static double kilo();
	static double hecto();
	static double deca();
	static double deci();
	static double centi();
	static double milli();
	static double micro();
	static double nano();
	static double pico();
	static double femto();
	static double atto();
	static double zepto();
	static double yocto();
};

#endif // CONVERTER_HPP_
