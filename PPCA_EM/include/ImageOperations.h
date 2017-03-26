//============================================================================
// Name        : ImageOperations.h
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================



#ifndef IMAGEOPERATIONS_H_
#define IMAGEOPERATIONS_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <tuple>
#include <math.h>
#include <sstream>
#include <armadillo>

using namespace arma;


void transform_Img(const double* inImg, double* outImg, int xpix, int ypix, double delta_x, double delta_y, double theta);
	/**Rotate and translate image.
	 * Arguments:
	 * 		inImg: const double*
	 * 			Input image.
	 * 		outImg: double*
	 * 			Output image. inImg and outImg are assumed to have the same dimensions and is not checked.
	 * 		xpix: int
	 * 			x pixels
	 * 		ypix: int
	 * 			y pixels
	 * 		delta_x: double
	 * 			Translational x shift in units of pixels.
	 * 		delta_y: double
	 * 			Translational y shift in units of pixels.
	 * 		theta: double
	 * 			Rotation shift in degrees.
	 * Returns:
	 *		void
	 */
void resample_Img(const double* inImg, double* outImg, int xpix, int ypix, int factor);
	/**	Resamples image by a given factor by zero padding and FFT.
	 * 	Arguments:
	 *		inImg: const double*
	 * 			Input image.
	 * 		outImg: double*
	 * 			Output image. inImg and outImg are assumed to have the same dimensions and is not checked.
	 * 		xpix: int
	 * 			x pixels
	 * 		ypix: int
	 * 			y pixels
	 * 		factor: int
	 * 			resample factor gas to be a factor of 2 or 1.
	 * 	Returns:
	 * 		void
	 * 	Throws:
	 * 		image_operations_exception
	 */

void rebin_Img(const double* inImg, double* outImg, int factor);
	/* Rebin image by a given factor by zero padding and FFT.
	 * Arguments:
	 *		inImg: const double*
	 * 			Input image.
	 * 		outImg: double*
	 * 			Output image. inImg and outImg are assumed to have the same dimensions and is not checked.
	 * 		xpix: int
	 * 			x pixels
	 * 		ypix: int
	 * 			y pixels
	 * 		factor: int
	 * 			rebin factor has to be a factor of 2 or 1.
	 * 	Returns:
	 * 		void
	 *	Throws:
	 * 		image_operations_exception
	 */





#endif /* IMAGEOPERATIONS_H_ */
