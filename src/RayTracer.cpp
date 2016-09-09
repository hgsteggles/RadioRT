/*
 * RayTracer.cpp
 *
 *  Created on: 12 Dec 2014
 *      Author: harry
 */

#include "RayTracer.hpp"
#include "Constants.hpp"
#include "Fluid.hpp"
#include "ProgressBar.hpp"
#include "RecombinationLine.hpp"

void RayTracerData::flip() {
	fluxFF.flip();
	tauFF.flip();
	emissionMeasureFF.flip();

	fluxRL.flip();
	tauRL.flip();
	emissionMeasureRL.flip();
}

TrigData::TrigData(double theta)
: sin_th(std::sin(theta))
, cos_th(std::cos(theta))
, tan_th(std::tan(theta))
, angle(theta)
{ }

double TrigData::sin() const { return sin_th; };
double TrigData::cos() const { return cos_th; };
double TrigData::tan() const { return tan_th; };
double TrigData::ang() const { return angle; };

void RayTracer::setImageDimensions(int xpixels, int ypixels) {
	pixels[0] = xpixels;
	pixels[1] = ypixels;
}

void RayTracer::setSampling(int s) {
	sampling = s;
}

void RayTracer::setFrequencies(const std::vector<double>& freqs) {
	frequencies = freqs;
}

void RayTracer::setIntegratingRL(bool on) {
	integratingRL = on;
}

void RayTracer::setIntegratingFF(bool on) {
	integratingFF = on;
}

void RayTracer::setDopplerShifted(bool dopplerShifted) {
	isDopplerShifted = dopplerShifted;
}

void RayTracer::setLineOfSightVelocity(double vLOS) {
	this->vLOS = vLOS;
}

void RayTracer::setDopplerIntMinPhiInc(double doppShiftPhiIncr) {
	this->doppShiftPhiIncr = doppShiftPhiIncr;
}

RayTracerData RayTracer::rayTrace3D(Fluid& fluid, double theta, double phi) {
	RayTracerData data;
	if (integratingFF) {
		data.fluxFF = DataCube(pixels[0], pixels[1], frequencies.size());
		data.tauFF = DataCube(pixels[0], pixels[1], frequencies.size());
		data.emissionMeasureFF = DataCube(pixels[0], pixels[1], frequencies.size());
	}
	if (integratingRL) {
		data.fluxRL = DataCube(pixels[0], pixels[1], frequencies.size());
		data.tauRL = DataCube(pixels[0], pixels[1], frequencies.size());
		data.emissionMeasureRL = DataCube(pixels[0], pixels[1], frequencies.size());
	}

	// Now perform the ray-tracing through the grid

	// If theta = 0.0, theta=1deg=pi/180.0
	// If theta = pi, theta=pi*(179.0/180.0)
	// (phi doesn't need these restrictions)

	if (theta == 0.0)
		theta = Constants::PI()/180.0; // 1deg in radians.
	if (theta == Constants::PI())
		theta = Constants::PI()*179.0/180.0; // 179deg in radians;

	// Calculate cos, sin and tan theta here, as well as cos and sin phi

	TrigData trig_th(theta);
	TrigData trig_ph(phi);

	// Calculate pixel size (given by fac) so that model grid fits into the image for any orientation. This ensures that the
	// image pixel size is independent of orientation.

	data.fac = 0.0;
	for (int i = 0; i < 3; ++i)
		data.fac += fluid.getNCells(i) * fluid.getNCells(i);
	data.fac = std::sqrt(data.fac) / (double)std::min(pixels[0], pixels[1]);

	double ds = 1.0 / sampling;

	int ist = fluid.getNCells(0) / 2;
	int jst = fluid.getNCells(1) / 2;
	int kst = fluid.getNCells(2) / 2;

	for (unsigned int ifreq = 0; ifreq < frequencies.size(); ++ifreq) {

		// Do the integration for an nu*nv array of lines of sight, uniformly spaced on the plane of the sky.
		// Loop through the horizontal pixels of the image.

		for (int iu = 0; iu < pixels[0]; ++iu) {
			// u coordinate of pixel in image plane.
			double u = data.fac * (double(iu - pixels[0] / 2) + 0.5);

			// loop through the vertical pixels of the image.
			for (int iv = 0; iv < pixels[1]; iv++){
				// v coordinate of pixel in image plane
				double v = data.fac * (double(iv - pixels[1] / 2) + 0.5);

				std::array<bool, 6> face;
				// Initialize faces logical
				for (int n = 0; n < 6; ++n)
					face[n] = false;

				// reset intensity, optical depth, emission measure.
				fluid.getBremsstrahlung().reset();
				fluid.getRecombinationLine().reset();

				// Determine the x,y,z coords of a position on the line of sight.
				// To do this we can imagine that the origin of the image plane
				// coincides with the origin of the hydro grid
				double ax = u * trig_ph.sin() - v * trig_th.cos() * trig_ph.cos();
				double ay = u * trig_ph.cos() + v * trig_th.cos() * trig_ph.sin();
				double az = v * trig_th.sin();

				// Now determine if the line of sight goes through the grid. Either it passes through 2 faces of the model grid, or through
				// no faces. First we set the components for the line of sight (x,y,z coords for a point on the line (ad,bd,cd), and x,y,z
				// components of the direction vector (add,bdd,cdd)).
				double ad = ax;               // )
				double bd = ay;               // )- Point on the los, ad = a' etc
				double cd = az;               // )
				double add = trig_ph.cos() * trig_th.sin();   //  )
				double bdd = trig_ph.sin() * trig_th.sin();   //  )- Direction vector of los, add = a" etc
				double cdd = trig_th.cos();          //  )

				std::array<double, 6> tp, xi, yi, zi;

				// For each face (plane), determine the intersection point with the line of sight

				double xmin = -(double)ist;
				double xmax = (double)(fluid.getNCells(0) - ist);
				double ymin = -(double)jst;
				double ymax = (double)(fluid.getNCells(1) - jst);
				double zmin = -(double)kst;
				double zmax = (double)(fluid.getNCells(2) - kst);

				std::array<double, 6> kp = std::array<double, 6>{xmax, xmin, ymax, ymin, zmax, zmin};

				tp[0] = (kp[0] - ad) / add;
				tp[1] = (kp[1] - ad) / add;
				tp[2] = (kp[2] - bd) / bdd;
				tp[3] = (kp[3] - bd) / bdd;
				tp[4] = (kp[4] - cd) / cdd;
				tp[5] = (kp[5] - cd) / cdd;
				for (int n = 0; n < 6; ++n) {
					xi[n] = ad + tp[n] * add;
					yi[n] = bd + tp[n] * bdd;
					zi[n] = cd + tp[n] * cdd;
				}

				// For each of these intersection points, determine whether the point
				// is within the bounds of each face. Modify the position by a small
				// amount such that the starting position is just within the grid bounds
				if ((yi[0] > ymin) && (yi[0] < ymax) && (zi[0] > zmin) && (zi[0] < zmax)){
					xi[0] = fluid.getNCells(0) - 1.0e-4;
					face[0] = true;
				}
				if ((yi[1] > ymin) && (yi[1] < ymax) && (zi[1] > zmin) && (zi[1] < zmax)){
					xi[1] = 0 + 1.0e-4;
					face[1] = true;
				}
				if ((xi[2] > xmin) && (xi[2] < xmax) && (zi[2] > zmin) && (zi[2] < zmax)){
					yi[2] = fluid.getNCells(1) - 1.0e-4;
					face[2] = true;
				}
				if ((xi[3] > xmin) && (xi[3] < xmax) && (zi[3] > xmin) && (zi[3] < zmax)){
					yi[3] = 0 + 1.0e-4;
					face[3] = true;
				}
				if ((xi[4] > xmin) && (xi[4] < xmax) && (yi[4] > ymin) && (yi[4] < ymax)){
					zi[4] = fluid.getNCells(2) - 1.0e-4;
					face[4] = true;
				}
				if ((xi[5] > xmin) && (xi[5] < xmax) && (yi[5] > ymin) && (yi[5] < ymax)){
					zi[5] = 0 + 1.0e-4;
					face[5] = true;
				}

				// Determine the number of intersections through the faces of the
				// model grid. There should be either 2 (los passes through grid)
				// or 0 (los misses grid).
				std::array<int, 6> facet;
				int ni = 0;
				for (int n = 0; n < 6; ++n){
					if (face[n]){
						facet[ni++] = n;
					}
				}

				if (ni == 0) {
					continue; // los misses model grid so next u,v pixel
				}
				if (ni != 2) {
					// NOTE: subsequent use of this code has revealed that on rare occasions 3 intersection points through the faces
					// are found. This is most likely due to numerical rounding. It might be possible to eliminate the incorrect value
					// in such cases.
					std::cout << "Unphysical number of intersections\n";
					std::cout << "iu = " << iu << ", iv = " << iv << "\n";
					std::cout << "ni = " << ni << "\n";
					for (int n = 0; n < ni; n++)
						std::cout << facet[n] << "\n";
					exit(0);
				}

				// Continue if there are two intersections. Determine at which intersection the line of sight enters the grid, and then
				// integrate from this point until the ray leaves the grid through the closer intersection point. Pick an imaginary point
				// on the los (suitably removed from the model grid) and determine the distance of the intersections from this. The
				// furthest intersection is where the los enters the model grid.
				double tim = 10.0 * std::max(fluid.getNCells(0), fluid.getNCells(1));
				double xim = ad + tim * add;
				double yim = bd + tim * bdd;
				double zim = cd + tim * cdd;
				double dist1 = std::pow(xi[facet[0]] - xim, 2) 
								+ std::pow(yi[facet[0]] - yim, 2) 
								+ std::pow(zi[facet[0]] - zim, 2);
				double dist2 = std::pow(xi[facet[1]] - xim, 2) 
								+ std::pow(yi[facet[1]] - yim, 2) 
								+ std::pow(zi[facet[1]] - zim, 2);

				// los enters through point 1 or point 2.
				double x = (dist2 < dist1) ? xi[facet[0]] : xi[facet[1]];
				double y = (dist2 < dist1) ? yi[facet[0]] : yi[facet[1]];
				double z = (dist2 < dist1) ? zi[facet[0]] : zi[facet[1]];

				// los wants to enter into the following grid cells (because we have normalized grid.ncells[0] etc by the cell size we do not need
				// to divide by the cell size here). Checked index is correct! (27/11/08)
				int ix = (x < 0.0) ? int(x) - 1 : int(x);
				int jy = (y < 0.0) ? int(y) - 1 : int(y);
				int kz = (z < 0.0) ? int(z) - 1 : int(z);

				int i = ix + ist;
				int j = jy + jst;
				int k = kz + kst;

				if ((i < 0) || (i >= fluid.getNCells(0)) || (j < 0) || (j >= fluid.getNCells(1)) || (k < 0) || (k >= fluid.getNCells(2))) {
					continue; // actually off grid
				}

				Cell const & cell = fluid.getCell(i, j, k);
				// now do the integration
				while (true) {
					// Determine emission and absorption coefficients for current cell
					if (integratingFF)
						fluid.getBremsstrahlung().update(cell, ds, frequencies[ifreq]);
					if (integratingRL)
						fluid.getRecombinationLine().update(cell, ds, frequencies[ifreq]);

					// calculate new position along line of sight, and new cell index
					x = x + ds * trig_ph.cos() * trig_th.sin();
					y = y + ds * trig_ph.sin() * trig_th.sin();
					z = z + ds * trig_th.cos();

					i = int(x) + ist;
					j = int(y) + jst;
					k = int(z) + kst;

					// has the l.o.s. reached the edge of the grid ? if so, set ok = false to terminate loop.

					if ((i < 0) || (i >= fluid.getNCells(0)) 
						|| (j < 0) || (j >= fluid.getNCells(1)) 
						|| (k < 0) || (k >= fluid.getNCells(2)))
						break;
				}

				// add the results of the last line of sight to the intensity and optical depth arrays

				if (integratingFF) {
					data.fluxFF(iu, iv, ifreq) = fluid.getBremsstrahlung().getIntensity();
					data.tauFF(iu, iv, ifreq) = fluid.getBremsstrahlung().getTau();
					data.emissionMeasureFF(iu, iv, ifreq) = fluid.getBremsstrahlung().getEmissionMeasure();
					data.intensityFF += fluid.getBremsstrahlung().getIntensity();
				}
				if (integratingRL) {
					data.fluxRL(iu, iv, ifreq) = fluid.getRecombinationLine().getIntensity();
					data.tauRL(iu, iv, ifreq) = fluid.getRecombinationLine().getTau();
					data.emissionMeasureRL(iu, iv, ifreq) = fluid.getRecombinationLine().getEmissionMeasure();
					data.intensityRL += fluid.getRecombinationLine().getIntensity();
				}
			} // v loop in image grid
		} // u loop in image grid
	} // freq loop in image grid

	return data;
}

std::array<int, 2> RayTracer::calcPixelNumberAxi(int rcells, int zcells, double dx, double desiredPixelSize) {
	std::array<int, 2> imageSize;

	double zeta = std::atan(0.5*zcells/rcells);
	// h is the largest possible extent in v of the projected grid for any theta between 0 and pi/2.
	// This is used to scale the grid so that the image always contains all of the grid, and the pixel size
	// will not vary if different thetas are used.
	double h = zcells*std::sin(zeta) + 2.0*rcells*std::cos(zeta);

	if (h/zcells < 1) {
		imageSize[0] = 2*1.01*rcells*dx/desiredPixelSize;
		imageSize[1] = (int)(imageSize[0]*zcells/(double)(rcells) + 0.5);
	}
	else {
		imageSize[1] = (int)(1.01*h*dx/desiredPixelSize + 0.5);
		imageSize[0] = (int)(imageSize[1]*rcells/(double)(zcells) + 0.5);
	}

	return imageSize;
}

double RayTracer::calcPixelSizeModificationAxi(int rcells, int zcells) {
	double zeta = std::atan(0.5*zcells/rcells);
	// h is the largest possible extent in v of the projected grid for any theta between 0 and pi/2.
	// This is used to scale the grid so that the image always contains all of the grid, and the pixel size
	// will not vary if different thetas are used.
	double h = zcells*std::sin(zeta) + 2.0*rcells*std::cos(zeta);
	int xpixels2 = pixels[0]/2;
	return 1.01*std::max(rcells/(double)xpixels2, h/(pixels[1]));
}

RayTracerData RayTracer::rayTraceAxiSymm(Fluid& fluid, double theta) {
	RayTracerData data;
	if (integratingFF) {
		data.fluxFF = DataCube(pixels[0], pixels[1], frequencies.size());
		data.tauFF = DataCube(pixels[0], pixels[1], frequencies.size());
		data.emissionMeasureFF = DataCube(pixels[0], pixels[1], frequencies.size());
	}
	if (integratingRL) {
		data.fluxRL = DataCube(pixels[0], pixels[1], frequencies.size());
		data.tauRL = DataCube(pixels[0], pixels[1], frequencies.size());
		data.emissionMeasureRL = DataCube(pixels[0], pixels[1], frequencies.size());
	}

	std::array<double, 2> nsamples;
	nsamples[0] = sampling*fluid.getNCells(0);
	nsamples[1] = sampling*fluid.getNCells(1);

	if (theta == 0.0)
		theta = Constants::PI()/180.0;
	if (theta == Constants::PI()/2.0)
		theta = Constants::PI()*89.0/180.0;

	TrigData trig_th(theta);

	int xpixels2 = pixels[0]/2;

	data.fac = RayTracer::calcPixelSizeModificationAxi(fluid.getNCells(0), fluid.getNCells(1));

	int jst = nsamples[1]*trig_th.cos()/(2.0*data.fac*sampling) - pixels[1]/2;

	ProgressBar progBar(frequencies.size()*xpixels2, 5, "RayTracing", false);

	for (unsigned int ifreq = 0; ifreq < frequencies.size(); ++ifreq) {
		for (int i = 1; i <= xpixels2; ++i) {
			double u = data.fac*sampling*(i - 0.5);
			if ( u < nsamples[0] ) {
				// Calculate limits of v for projected model grid.
				double sin_thsqrt = trig_th.sin()*std::sqrt((nsamples[0] - u)*(nsamples[0] + u));
				double vtop = nsamples[1]*trig_th.cos() + sin_thsqrt;
				double vbot = -sin_thsqrt;

				// Loop through vertical pixels.
				for (int j = 1; j <= pixels[1]; ++j) {
					double v = data.fac*sampling*(j + jst - 0.5);

					// Continue if LOS goes through the grid.
					if (v < vtop && v > vbot) {
						fluid.getBremsstrahlung().reset();
						fluid.getRecombinationLine().reset();

						// calculate point of ingress on grid
						// phir is the point at which LOS crosses cylinder r = R
						// phiz is the point at which LOS crosses plane z = 0
						double phir = Constants::PI() - std::asin(u/nsamples[0]);
						double phiz = - std::atan(u*trig_th.sin()/v);

						// Make sure phiz is in the right quadrant.
						if (phiz < 0.0)
							phiz += Constants::PI();

						double phi = std::min(phir, phiz);
						double r = u/std::sin(phi);
						double z = z = v/trig_th.cos() + u*trig_th.tan()/std::tan(phi);

						int ic, jc;
						if (phiz < phir) {
							// LOS enters model grid through plane z = 0.
							ic = int(r/sampling) + 1;
							jc = 1;
						}
						else {
							phi = phir;
							ic = fluid.getNCells(0);
							jc = int(z/sampling) + 1;
						}

						// Now do integration.
						while (true) {
							int cellx = ic-1;
							int celly = jc-1;

							// Calculate value of phi corresponding to the point at which the LOS enters the next grid cell in r and z.
							// inr is the next cell in r (direction depends on whether phi>pi/2 or pi/2.
							// inz is the next cell in z (always up given theta > 0).
							int inr;
							if (phi > Constants::PI()/2.0) {
								inr = ic - 1;

								// Check for the special case where the LOS crosses the line phi = pi/2 inside the current cell. In this case the LOS
								// crosses the same cell boundary twice and the next cell boundary in r is the same as the last cell boundary in r.

								if (inr*sampling <= u) { // LOS crosses phi=pi/2 in the current cell.
									inr = ic;
									phir = std::asin(u/(sampling*inr));
								}
								else { // LOS doesn't cross phi=pi/2 in this cell.
									phir = Constants::PI() - std::asin(u/(sampling*inr));
								}
							}
							else {
								inr = ic + 1;
								phir = std::asin(u/(sampling*inr));
							}

							int inz = jc + 1;
							phiz = std::atan(u*trig_th.tan()/((inz-1)*sampling - v/trig_th.cos()));

							// Make sure phiz is in the right quadrant.
							if (phiz < 0.0)
								phiz += Constants::PI();

							double prev_phi = phi;
							// Whichever of phir, phiz has the smallest value decides whether the LOS next crosses a z-boundary or an r-boundary.
							if (phir < phiz) {
								++jc; // Move to next cell in z, stay in same r cell.
								phi = phiz; // Set current phi to new phi.
							}
							else {
								ic = inr; // Move to next cell in r, stay in same z cell.
								phi = phir; // Set current phi to new phi.
							}

							double ds = std::abs((u/trig_th.cos())*(1.0/std::tan(phi) - 1.0/std::tan(prev_phi))); // Calculate path-length through last cell.

							Cell const & cell = fluid.getCell(cellx, celly, 0);

							if (isDopplerShifted) {
								double ur = cell.velocity[0];
								double uz = cell.velocity[1];
								DopplerShift dopp(frequencies[ifreq], vLOS, ur, uz, u, trig_th, prev_phi, phi, ds, doppShiftPhiIncr);

								for (int idopp = 0; idopp < dopp.deltas.size(); ++idopp) {
									if (integratingFF)
										fluid.getBremsstrahlung().update(cell, dopp.deltas[idopp], frequencies[ifreq]*dopp.shifts[idopp]);
									if (integratingRL)
										fluid.getRecombinationLine().update(cell, dopp.deltas[idopp], frequencies[ifreq]*dopp.shifts[idopp]);

								}
							}
							else {
								if (integratingFF)
									fluid.getBremsstrahlung().update(cell, ds, frequencies[ifreq]);
								if (integratingRL)
									fluid.getRecombinationLine().update(cell, ds, frequencies[ifreq]);

							}

							// Calculate r and z corresponding to new cell-boundary crossing point.
							r = u / std::sin(phi);
							z = v / trig_th.cos() + u * trig_th.tan() / std::tan(phi);

							// LOS reached the edge of the grid ? then ok = false.
							if ((ic >= fluid.getNCells(0) && phi < Constants::PI()/2.0) || jc >= fluid.getNCells(1))
								break;
						}

						// Add the results of the last LOS to the intensity and optical depth arrays.
						if (integratingFF) {
							data.fluxFF(i + xpixels2 - 1, j - 1, ifreq) = fluid.getBremsstrahlung().getIntensity();
							data.fluxFF(xpixels2 - i, j - 1, ifreq) = fluid.getBremsstrahlung().getIntensity();

							data.tauFF(i + xpixels2 - 1, j - 1, ifreq) = fluid.getBremsstrahlung().getTau();
							data.tauFF(xpixels2 - i, j - 1, ifreq) = fluid.getBremsstrahlung().getTau();

							data.emissionMeasureFF(i + xpixels2 - 1, j - 1, ifreq) = fluid.getBremsstrahlung().getEmissionMeasure();
							data.emissionMeasureFF(xpixels2 - i, j - 1, ifreq) = fluid.getBremsstrahlung().getEmissionMeasure();

							data.intensityFF += 2.0*fluid.getBremsstrahlung().getIntensity();
						}

						if (integratingRL) {
							data.fluxRL(i + xpixels2 - 1, j - 1, ifreq) = fluid.getRecombinationLine().getIntensity();
							data.fluxRL(xpixels2 - i, j - 1, ifreq) = fluid.getRecombinationLine().getIntensity();

							data.tauRL(i + xpixels2 - 1, j - 1, ifreq) = fluid.getRecombinationLine().getTau();
							data.tauRL(xpixels2 - i, j - 1, ifreq) = fluid.getRecombinationLine().getTau();

							data.emissionMeasureRL(i + xpixels2 - 1, j - 1, ifreq) = fluid.getRecombinationLine().getEmissionMeasure();
							data.emissionMeasureRL(xpixels2 - i, j - 1, ifreq) = fluid.getRecombinationLine().getEmissionMeasure();

							data.intensityRL += 2.0 * fluid.getRecombinationLine().getIntensity();
						}
					}
				}
			}
			progBar.update(ifreq*xpixels2 + i);
		}
	}

	return data;
}

DopplerShift::DopplerShift(double frequency, double vLOS, double vr, double vz, double x, const TrigData& trig_th, double phi1, double phi2, double path_length, double min_angle)
{
		double c = Constants::lightSpeed();

		double u = std::sqrt(vz*vz + vr*vr);

		double phi = phi1;
		double ds = 0;

		while (phi > phi2) {
			double cos_ph = std::cos(phi);
			double cosa = (vr*cos_ph + vz*trig_th.cos())/(u*std::sqrt(cos_ph*cos_ph + trig_th.cos()*trig_th.cos()));

			double min_phi = std::min(min_angle, std::abs(phi2 - phi));
			//double min_ds = radius*min_phi/(std::sin(phi)*trig_th.sin());
			double min_ds = std::abs((x/trig_th.cos())*(1.0/std::tan(phi-min_phi) - 1.0/std::tan(phi)));
			ds += min_ds;
			deltas.push_back(min_ds);

			double shift = std::sqrt((1.0 - (u*cosa - vLOS)/c)/(1.0 + (u*cosa - vLOS)/c));

			shifts.push_back(shift);

			phi -= min_phi;
		}
	}





