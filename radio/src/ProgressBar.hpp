/** Provides the ProgressBar class.
 *
 * @file ProgressBar.hpp
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014 - the first version.
 */

#ifndef PROGRESSBAR_HPP_
#define PROGRESSBAR_HPP_

#include "Timer.hpp"

#include <stdint.h>
#include <chrono>
#include <string>

static double dummy_checkpoint;

/**
 * @class ProgressBar
 *
 * @brief UI feature: an indicator of progress is displayed in the terminal.
 *
 * @version 0.8, 24/11/2014
 */
class ProgressBar
{
public:
	static int messageSpace;
	ProgressBar(double tmax, int cpoint, const std::string msg, bool debug);
	bool update(double timeCurrent, bool output_on = true, double& dt_nextCheckpoint = dummy_checkpoint);
	void end(bool output_on);
	void reset(double tmax, int cpoint, const std::string msg);
private:
	Timer timer;
	typedef std::chrono::steady_clock Clock;
	typedef Clock::time_point time_point;
	typedef Clock::period period;
	typedef std::chrono::duration<float, period> duration;
	double timeTotal;
	int checkpoint;
	std::string messageProgress;
	double percentProgress;
	int checkpointProgress;
	time_point clockLast;
	duration speedEMA;
	double smoothing;
	bool isStarting;
	bool debugging;
};

#endif //PROGRESSBAR_HPP_
