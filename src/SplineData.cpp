#include "SplineData.hpp"

#include <cmath>
#include <stdexcept>
#include <iostream>

SplineData2D::SplineData2D()
: m_empty(true)
{

}

SplineData2D::SplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f)
: m_x(x)
, m_y(y)
, m_f(f)
, m_fftmp(x.size(), 0)
, m_empty(false)
{
	if (f.size() < 1 || f.size() != x.size() || f[0].size() != y.size())
		throw std::runtime_error("SplineData2D::SplineData2D: independent variable vector size incompatible with tabulated function grid size.");
}

bool SplineData2D::isEmpty() {
	return m_empty;
}

/**
 * @brief Calculates second derivatives of an interpolating function.
 * @param x Vector of function arguments.
 * @param f Vector of function values given a vector of function arguments.
 * @param yp1 First derivative of function at x[0].
 * @param ypn First derivative of function at x[n-1].
 * @return Vector of second derivatives of interpolating function.
 */
std::vector<double> SplineData2D::spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn) {
	int n = x.size();
	std::vector<double> y2(n, 0.0);
	std::vector<double> u(n-1, 0.0);
	int i,k;
	double p,qn,sig,un;
	//The lower boundary condition is set either to be "natural"
	//or else to have a specified first derivative.
	if ( yp1 > 0.99e30) {
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((f[1]-f[0])/(x[1]-x[0])-yp1);
	}
	//This is the decomposition loop of the tridiagonal algorithm.
	//y2 and u are used for temporary storage of the decomposed factors.
	for ( i = 1; i < n-1; i++ ) {
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (f[i+1]-f[i])/(x[i+1]-x[i]) - (f[i]-f[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	//The upper boundary condition is set either to be "natural"
	// or else to have a specified first derivative.
	if ( ypn > 0.99e30 ) {
		qn = 0.0;
		un = 0.0;
	}
	else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(f[n-1]-f[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	//This is the backsubstitution loop of the tridiagonal algorithm.
	for ( k = n-2; k >= 0; k-- )
		y2[k] = y2[k]*y2[k+1]+u[k];
	return y2;
}

double SplineData2D::splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2) {
	int klo, khi, k;
	double h,b,a;
	int n = x.size();
	//We will find the right place in the table by means of bisection. This is
	//optimal if sequential calls to this routine are at random values of x.
	//If sequential calls are in order, and closely spaced, one would do better
	//to store previous values of klo and khi and test if they remain appropriate
	//on the next call.
	klo = 0;
	khi = n-1;
	while(khi-klo > 1){
		k = (khi+klo) >> 1;
		if(x[k] > x2)
			khi = k;
		else
			klo = k;
	}
	//klo and khi now bracket the input value of x.
	h = x[khi]-x[klo];
	a = (x[khi]-x2)/h;
	b = (x2-x[klo])/h; //Cubic spline polynomial is now evaluated.

	return a*f[klo]+b*f[khi]+((a*a*a-a)*f2[klo]+(b*b*b-b)*f2[khi])*(h*h)/6.0;
}

std::vector<std::vector<double>> SplineData2D::splie2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f) {
	int m = x.size();
	int n = y.size();

	if (m == 0)
		throw std::runtime_error("SpineData2D::splie2: x vector has 0 elements");
	if (n == 0)
		throw std::runtime_error("SpineData2D::splie2: y vector has 0 elements");

	std::vector<std::vector<double>> f2;

	for (int j = 0; j < m; ++j)
		f2.push_back(spline(y, f[j], 1.0e30, 1.0e30));

	return f2;
}

double SplineData2D::splin2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f,
		const std::vector<std::vector<double>>& f2, double x2, double y2) {
	int m = x.size();
	int n = y.size();

	for (int j = 0; j < m; ++j)
		m_fftmp[j] = splint(y, f[j], f2[j], y2);

	return splint(x, m_fftmp, spline(x, m_fftmp, 1.0e30, 1.0e30), x2);
}

LinearSplineData2D::LinearSplineData2D()
: SplineData2D()
{

}

LinearSplineData2D::LinearSplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f)
: SplineData2D(x, y, f)
{
	m_f2 = splie2(m_x, m_y, m_f);
}

/**
 * @brief Calculates a cubic spline interpolated value from the data this object was initialised with.
 * @param x Interpolation location.
 * @return Cubic spline interpolated value.
 */
double LinearSplineData2D::interpolate(double x, double y) {
	double rate = 0;
	if (x > m_x[m_x.size()-1] || x < m_x[0] || y > m_y[m_y.size()-1] || y < m_y[0])
		throw std::runtime_error("LinearSplineData2D::interpolate: out of bounds (x = " + std::to_string(x) + ", y = " + std::to_string(y) + ", xminmax = [" + std::to_string(m_x[0]) + ", " + std::to_string(m_x[m_x.size()-1]) + "], yminmax = [" + std::to_string(m_y[0]) + ", " + std::to_string(m_y[m_y.size()-1]) + "]");
	else
		rate = splin2(m_x, m_y, m_f, m_f2, x, y);
	return rate;
}

LogSplineData2D::LogSplineData2D()
: SplineData2D()
{

}

LogSplineData2D::LogSplineData2D(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f)
: SplineData2D(x, y, f)
{
	m_f2 = splie2(m_x, m_y, m_f);
}

double LogSplineData2D::interpolate(double x, double y) {
	double rate = 0;
	if (x > m_x[m_x.size()-1] || x < m_x[0] || y > m_y[m_y.size()-1] || y < m_y[0])
		throw std::runtime_error("LinearSplineData2D::interpolate: out of bounds (x = " + std::to_string(x) + ", y = " + std::to_string(y) + ", xminmax = [" + std::to_string(m_x[0]) + ", " + std::to_string(m_x[m_x.size()-1]) + "], yminmax = [" + std::to_string(m_y[0]) + ", " + std::to_string(m_y[m_y.size()-1]) + "]");
	else
		rate = splin2(m_x, m_y, m_f, m_f2, x, y);
	return rate;
}
