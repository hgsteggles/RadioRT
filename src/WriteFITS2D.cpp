/*
 * WriteFITS2D.cpp
 *
 *  Created on: 18 Dec 2014
 *      Author: harry
 */

//
// subroutine WFITS2D
//
// Routine for writing out a 2-D array of data z, of dimensions
// nx by ny into a FITS format file. The output filename is
// filename.fits.
//
// Author: J.M. Pittard
//
// Original version: 27/11/08 - based on previous radio and X-ray codes

#undef MAIN

#include "WriteFITS2D.hpp"

#include "cfitsio/fitsio.h"      //  For FITSIO (from cfitsio directory)
#include "cfitsio/longnam.h"     //  For FITSIO (from cfitsio directory)

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

void WriteFITS::wfits(std::string filename, std::string units, const DataCube& cube, RadioRT_Parameters params,
		   double pixdeg, double pixsize, double dx, double total) {
	double *z = new double[cube.sizex()*cube.sizey()*cube.sizez()];

	double image_max = 0;
	double image_min = cube(0, 0, 0);
	for (unsigned int k = 0, num = 0; k < cube.sizez(); ++k) {
		for (unsigned int j = 0; j < cube.sizey(); ++j) {
			for (unsigned int i = 0; i < cube.sizex(); ++i) {
				z[num++] = cube(i, j, k);

				image_max = std::max(image_max, cube(i, j, k));
				image_min = std::min(image_min, cube(i, j, k));
			}
		}
	}

	int naxis = cube.sizez() == 1 ? 2 : 3;

	FILE *fptr;
	fitsfile *fitsptr;                      //  Pointer to fits file.
	int status = 0;                         //  FITSIO status.
	int bitpix = FLOAT_IMG;                 //  Single precision.

	filename = filename + ".fits";

	remove( filename.c_str() );
	if ( fits_create_file(&fitsptr, filename.c_str(), &status) ){
		WriteFITS::printFITS_Error( status );
	}

	//Write required header info
	long naxes[3];
	naxes[0] = cube.sizex();
	naxes[1] = cube.sizey();
	naxes[2] = cube.sizez();

	if ( fits_create_img(fitsptr, bitpix, naxis, naxes, &status) ){
		std::cerr << "  Error creating header info.\n";
		WriteFITS::printFITS_Error( status );
	}

	if ( fits_write_date(fitsptr, &status ) ){
		std::cerr << "  Error writing DATE keyword.\n";
		WriteFITS::printFITS_Error( status );
	}

	char ctype1[FLEN_VALUE];  //  X-axis label
	char ctype2[FLEN_VALUE];  //  Y-axis label
	char ctype3[FLEN_VALUE];  //  Z-axis label
	char ctype4[FLEN_VALUE];  //  S-axis label
	strcpy(ctype1,"RA---SIN");
	strcpy(ctype2,"DEC--SIN");
	strcpy(ctype3,"FREQ");
	strcpy(ctype4,"STOKES");

	double pixs = pixdeg;
	double mpixs = -pixdeg;
	double crpix1 = double(cube.sizex()/2);
	double crpix2 = double(cube.sizey()/2);

	double channelWidth = params.bandwidth/params.nchannels;

	fits_update_key(fitsptr, TSTRING, "CTYPE1", ctype1, "Axis type", &status);
	fits_update_key(fitsptr, TDOUBLE, "CDELT1", &mpixs, "Axis coordinate increment (deg)", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX1", &crpix1, "Axis coordinate reference pixel",  &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL1", &params.rightAscension, "Axis coordinate value at CRPIX", &status);
	fits_update_key(fitsptr, TSTRING, "CTYPE2", ctype2, "Axis type", &status);
	fits_update_key(fitsptr, TDOUBLE, "CDELT2", &pixs, "Axis coordinate increment (deg)", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX2", &crpix2, "Axis coordinate reference pixel", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL2", &params.declination, "Axis coordinate value at CRPIX", &status);

	double channelRef = params.nchannels / 2.0;

	fits_update_key(fitsptr, TSTRING, "CTYPE3", ctype3,"Axis type", &status);
	fits_update_key(fitsptr, TDOUBLE, "CDELT3", &channelWidth,"Channel width (Hz)", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX3", &channelRef,"Axi coordinate reference channel", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL3", &params.frequency,"Obs Frequency (Hz)", &status);

	double stokes = params.stokes;

	fits_update_key(fitsptr, TSTRING, "CTYPE4", ctype4, "Axis type", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRVAL4", &stokes, "Unity", &status);
	fits_update_key(fitsptr, TDOUBLE, "CRPIX4", &stokes, "Unity", &status);
	fits_update_key(fitsptr, TDOUBLE, "CDELT4", &stokes, "Unity", &status);

	double epoch = 2000.000;
	double bpa  = 0.0;
	fits_update_key(fitsptr, TDOUBLE, "EPOCH", &epoch, "EPOCH", &status);
	fits_update_key(fitsptr, TDOUBLE, "EQUINOX", &epoch, "EQUINOX", &status);
	fits_update_key(fitsptr, TDOUBLE, "BPA", &bpa, "PA (radians)", &status);

	char * units_writable = new char[units.size() + 1];
	std::copy(units.begin(), units.end(), units_writable);
	units_writable[units.size()] = '\0'; // don't forget the terminating 0
	fits_update_key(fitsptr,TSTRING,"BUNIT", units_writable, "Units", &status);

	fits_update_key(fitsptr, TINT, "Samples", &params.sampling, "Samples per cell", &status);
	fits_update_key(fitsptr, TDOUBLE, "Theta", &params.theta, "Inclination (deg)", &status);
	fits_update_key(fitsptr, TDOUBLE, "Phi", &params.phi, "Azimuth (deg)", &status);
	fits_update_key(fitsptr, TDOUBLE, "RESTFREQ", &params.frequency, "Frequency (Hz)", &status);
	fits_update_key(fitsptr, TDOUBLE, "Dist", &params.dist, "Distance (kpc)", &status);
	fits_update_key(fitsptr, TDOUBLE, "PIXCM", &pixsize, "Pixel size (cm)", &status);
	double pixarc = pixdeg * 60.0 * 60.0;
	fits_update_key(fitsptr, TDOUBLE, "PIXAS", &pixarc, "Axes coordinate increment (arcsec)", &status);
	std::string mpix_label = "Brightest Pixel (" + units + ")";
	fits_update_key(fitsptr, TDOUBLE, "DATAMIN", &image_min, mpix_label.c_str(), &status);
	fits_update_key(fitsptr, TDOUBLE, "DATAMAX", &image_max, mpix_label.c_str(), &status);
	if (total >= 0) {
		fits_update_key(fitsptr, TDOUBLE, "TOTAL", &total, "Total", &status);
	}

	//  Write the 1-D array, cunningly stored in the same way as a 2-D
	//  image, into the FITS image.
	long group = 1;
	if (naxis == 2) {
		if ( fits_write_2d_dbl(fitsptr, group, naxes[0], naxes[0], naxes[1], &z[0], &status) )
			WriteFITS::printFITS_Error( status );
	}
	else {
		if ( fits_write_3d_dbl(fitsptr, group, naxes[0], naxes[1], naxes[0], naxes[1], naxes[2], &z[0], &status) )
			WriteFITS::printFITS_Error( status );
	}

	//  Close the output fits file.
	if ( fits_close_file(fitsptr, &status) )  //  Close the file.
		WriteFITS::printFITS_Error( status );

	delete []z;
}

//file FitsPrintError.c
/*
 + Print out cfitsio error messages and exit program.
   Taken directly from FITSIO cookbook.c
 */

/**
 * @brief Prints cfitsio error messages and exits program.
 * @param status Error code.
 */
void WriteFITS::printFITS_Error( int status ){
	if (status) {
		fits_report_error(stderr, status);    //  Print error report.
		throw std::runtime_error("FITS error code = " + std::to_string(status));
	}
}






