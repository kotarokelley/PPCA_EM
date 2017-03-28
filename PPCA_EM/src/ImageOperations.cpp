//============================================================================
// Name        : ImageOperations.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include "ImageOperations.h"

void transform_Img_single(const double* inImg, double* outImg, int xpix, int ypix, double delta_x, double delta_y, double theta, double max_radius){

	double radians = theta * M_PI / 180.0;

	mat trans_mat(3,3,fill::zeros);

	trans_mat(0,0) = std::cos(radians);
	trans_mat(0,1) = std::sin(radians);
	trans_mat(0,2) = delta_x;
	trans_mat(1,0) = -std::sin(radians);
	trans_mat(1,1) = std::cos(radians);
	trans_mat(1,2) = delta_y;
	trans_mat(2,2) = 1.0;

	for (int i=0; i<xpix; i++){

		for (int j=0; j<ypix; j++){

			vec coords_old(3);
			coords_old(0) = i-xpix/2;					// Set center of image at xpix/2 and ypix/2
			coords_old(1) = j-ypix/2;
			coords_old(2) = 1;

			vec coords_new = trans_mat * coords_old;	// Apply transformation matrix.

			double x_new = coords_new(0) + xpix/2;		// Re set center image at origin (0,0).
			double y_new = coords_new(1) + ypix/2;

			double mean(0);

			for (int m=0; m<xpix*ypix; m++)
				mean += inImg[m];

			mean /= (xpix*ypix);

			if (x_new>=xpix || y_new>=ypix || std::sqrt(std::pow(x_new/2,2) + std::pow(y_new/2,2)) > max_radius)
				outImg[i*xpix + j] = mean;

			else{
				double x2x1, y2y1, x2x, y2y, yy1, xx1;

				int x1 = std::floor(x_new), x2 = std::ceil(x_new), y1 = std::ceil(y_new), y2 = std::floor(y_new);

				if (x1== x2)						// This takes care of landing exactly on a grid point.
					x2 +=1;
				if (y2 == y1)
					y1 +=1;

				x2x1 = x2 - x1;
				y2y1 = y2 - y1;
				x2x = x2 - x_new;
				y2y = y2 - y_new;
				yy1 = y_new - y1;
				xx1 = x_new - x1;

				outImg[j*xpix + i] = //inImg[y1*xpix + x1];
					1.0 / (x2x1 * y2y1) *
					(inImg[y1*xpix + x1] * x2x * y2y +
					inImg[y1*xpix + x2] * xx1 * y2y +
					inImg[y2*xpix + x1] * x2x * yy1 +
					inImg[y2*xpix + x2] * xx1 * yy1);

			}
		}
	}
}

void transform_Img_multiple(const double* inImg, double* outImg, int xpix, int ypix, int numImg, double delta_x, double delta_y, double theta, double max_radius){

	mat trans_mat(3,3,fill::zeros);
	trans_mat(0,0) = std::cos(theta);
	trans_mat(0,1) = -std::sin(theta);
	trans_mat(0,2) = delta_x;
	trans_mat(1,0) = std::sin(theta);
	trans_mat(1,1) = std::cos(theta);
	trans_mat(1,2) = delta_y;
	trans_mat(2,2) = 1.0;


	double* mean = new double[numImg]();					// Pre-compute the mean of each image to be used for out of bounds conditions.
	for (int k=0; k<numImg; k++){
		int offset = k*xpix*ypix;
		for (int m=0; m<xpix*ypix; m++)
			mean[k] += inImg[m+offset];
		mean[k] /= xpix*ypix;
	}

	for (int i=0; i<xpix; i++){
		for (int j=0; j<ypix; j++){
			vec coords_old(3);
			coords_old(0) = i-xpix/2;					// Set center of image at xpix/2 and ypix/2
			coords_old(1) = j-ypix/2;
			coords_old(2) = 1;

			vec coords_new = trans_mat * coords_old;	// Apply transformation matrix.

			double x_new = coords_new(0) + xpix/2;		// Re set center image at origin (0,0).
			double y_new = coords_new(0) + ypix/2;

			if (x_new>=xpix || y_new>=ypix || std::sqrt(std::pow(x_new/2,2) + std::pow(y_new/2,2)) > max_radius){
				for (int k=1; k < numImg; k++)
					outImg[i*xpix + k*xpix*ypix] = mean[k];
			}

			else{
				double x2x1, y2y1, x2x, y2y, yy1, xx1;
				int x1 = std::floor(x_new), x2 = std::ceil(x_new), y1 = std::floor(y_new), y2 = std::ceil(y_new);
				x2x1 = x2 - x1;
				y2y1 = y2 - y1;
				x2x = x2 - x_new;
				y2y = y2 - y_new;
				yy1 = y_new - y1;
				xx1 = x_new - x1;
				for (int k=1; k < numImg; k++){
					int offset = k*xpix*ypix;
					outImg[i*xpix + j + offset] =
						1.0 / (x2x1 * y2y1) *
						(inImg[x1*xpix + y1 + offset] * x2x * y2y +
						inImg[x2*xpix + y1 + offset] * xx1 * y2y +
						inImg[x1*xpix + y2 + offset] * x2x * yy1 +
						inImg[x2*xpix + y2 + offset] * xx1 * yy1);
				}
			}
		}
	}

}
void resample_Img(const double* inImg, double* outImg, int xpix, int ypix, int factor){

}


void rebin_Img(const double* inImg, double* outImg, int factor){

}
