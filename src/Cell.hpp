//
// Created by harry on 17/02/16.
//

#ifndef RADIO_CELL_HPP
#define RADIO_CELL_HPP

struct SPECIES {
	enum ID { H, HE, CNO, N };
};

class Cell {
public:
	Cell() {
		for (int i = 0; i < SPECIES::N; ++i)
			massFractions[i] = 0;
	}

	double dx = 0;

	double density = 0;
	double pressure = 0;
	double hii = 0;
	std::array<double, 3> velocity = std::array<double, 3>{0, 0, 0};
	double temperature = 0;
	double massFractionH = 1.0;

	std::array<double, SPECIES::N> massFractions;

	double absorption_ff_coeff = 0;
	double emission_ff_coeff = 0;
	double absorption_l_coeff = 0;
	double sigma = 0;

	double bn = 0;
	double bnp1 = 0;
};

#endif //RADIO_CELL_HPP
