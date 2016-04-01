/** Provides the SplineData class.
 *
 * @file SplineData.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef SPLINEDATA_HPP_
#define SPLINEDATA_HPP_

#include <utility>
#include <vector>

/**
 * @class SplineData2D
 *
 * @brief Contains cubic spline data for interpolation.
 *
 * Cubic spline interpolation and powerlaw extrapolation (using logarithm gradients at the extremes).
 *
 * @version 0.8, 24/11/2014
 */
class SplineData2D {
public:
	SplineData2D();
	SplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
	virtual ~SplineData2D() {};

	bool isEmpty();

	virtual double interpolate(double x, double y) = 0;
	static std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn);
	static double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2);
	std::vector<std::vector<double>> splie2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
	double splin2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f,
			const std::vector<std::vector<double>>& f2, double x2, double y2);
protected:
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_fftmp;
	std::vector<std::vector<double>> m_f;
	std::vector<std::vector<double>> m_f2;

	bool m_empty;
};

class LinearSplineData2D : public SplineData2D {
public:
	LinearSplineData2D();
	LinearSplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
	~LinearSplineData2D() {};

	virtual double interpolate(double x, double y);
};

class LogSplineData2D : public SplineData2D {
public:
	LogSplineData2D();
	LogSplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
	~LogSplineData2D() {};

	virtual double interpolate(double x, double y);
};



#endif // SPLINEDATA_HPP_
