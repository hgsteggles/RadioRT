/*
 * RadioRT.hpp
 *
 *  Created on: 8 Jan 2015
 *      Author: harry
 */

#ifndef RADIORT_HPP_
#define RADIORT_HPP_

#include "RadioRT_Parameters.hpp"

class RadioRT {
public:
	RadioRT(const RadioRT_Parameters& p);
	void run();
private:
	RadioRT_Parameters params;
};



#endif /* RADIORT_HPP_ */
