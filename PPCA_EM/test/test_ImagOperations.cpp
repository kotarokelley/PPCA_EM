
#include "ImageOperations.h"
#include "Parser.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui.hpp"

int main(void){

	char * filename = "zip_test.mrc";

	mrcParser testParser = mrcParser(filename);  // input binary data

	char * f_data = testParser.getData(false);			// extract data, verbose off

	char * f_header = testParser.getHdr();

	std::vector<int> f_dim = testParser.getDim();					// get dimension of image stack TODO, this needs to be fixed since when testParser goes out of scope, f_dim will go out of scope.

	int f_pixelType = testParser.getPixelType();		// get pixel type


	std::cout << "Data Set:\n";
	std::cout << filename << "\n";
	std::cout << "pixel type\n";
	std::cout << f_pixelType << "\n";

	std::cout << "dimensions of data set\n";
	std::cout << f_dim[0] << "\n";
	std::cout << f_dim[1] << "\n";
	std::cout << f_dim[2] << "\n";

	std::cout << "Casting raw data array to float array.\n";

	float * f_data_float = new float[f_dim[0]*f_dim[1]*f_dim[2]];

	std::memcpy(f_data_float, f_data, sizeof(float)*f_dim[0]*f_dim[1]*f_dim[2]);

	double * f_data_double = new double[f_dim[0]*f_dim[1]*f_dim[2]];
	std::cout << "\nConverting data to double for input into pca object.\n";
	for (int i=0; i<f_dim[0]*f_dim[1]*f_dim[2]; i++){
		f_data_double[i] = f_data_float[i];
	}

	std::cout << "Testing function transform_Img_single\n";
	std::cout << "Displaying first image\n";
	cv::Mat img_0(cv::Size(f_dim[0],f_dim[1]),CV_64F,f_data_double);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_0);
	cv::waitKey(2000);

	double * rotated_img_0 = new double[f_dim[0]*f_dim[1]]();
	transform_Img_single(f_data_double, rotated_img_0, f_dim[0], f_dim[1], 0,  0,  90,  f_dim[0]/2 + 4);

	std::cout << "Displaying first image after transformation\n";
	cv::Mat img_0_rot(cv::Size(f_dim[0],f_dim[1]),CV_64F,rotated_img_0);
	cv::imshow("image_0_rot", img_0_rot);
	cv::waitKey(2000);

	std::cout << "\nTesting function transform_Img_multiple\n";
	std::cout << "Displaying a few images\n";
	cv::Mat img_5(cv::Size(f_dim[0],f_dim[1]),CV_64F, &f_data_double[f_dim[0]*f_dim[1]*5]);
	cv::Mat img_200(cv::Size(f_dim[0],f_dim[1]),CV_64F,&f_data_double[f_dim[0]*f_dim[1]*200]);
	cv::Mat img_1000(cv::Size(f_dim[0],f_dim[1]),CV_64F,&f_data_double[f_dim[0]*f_dim[1]*1000]);
	cv::Mat win_mat(cv::Size(f_dim[0]*3,f_dim[1]), CV_64F);
	img_5.copyTo(win_mat(cv::Rect(  0, 0, f_dim[0], f_dim[1])));
	img_200.copyTo(win_mat(cv::Rect( f_dim[0], 0, f_dim[0], f_dim[1])));
	img_1000.copyTo(win_mat(cv::Rect( f_dim[0]*2, 0, f_dim[0], f_dim[1])));

    cv::namedWindow("Images");
	cv::imshow("Images", win_mat);
	cv::waitKey(3000);

	std::cout << "Displaying a few images after transformation\n";

	double * imgs_rotated = new double[f_dim[0]*f_dim[1]*f_dim[2]]();
	transform_Img_multiple(f_data_double, imgs_rotated, f_dim[0], f_dim[1], f_dim[2], 0,  0,  90,  f_dim[0]/2 + 4);
	cv::Mat img_5_rot(cv::Size(f_dim[0],f_dim[1]),CV_64F, &imgs_rotated[f_dim[0]*f_dim[1]*5]);
	cv::Mat img_200_rot(cv::Size(f_dim[0],f_dim[1]),CV_64F,&imgs_rotated[f_dim[0]*f_dim[1]*200]);
	cv::Mat img_1000_rot(cv::Size(f_dim[0],f_dim[1]),CV_64F,&imgs_rotated[f_dim[0]*f_dim[1]*1000]);
	cv::Mat win_mat_rot(cv::Size(f_dim[0]*3,f_dim[1]), CV_64F);
	img_5_rot.copyTo(win_mat_rot(cv::Rect(  0, 0, f_dim[0], f_dim[1])));
	img_200_rot.copyTo(win_mat_rot(cv::Rect( f_dim[0], 0, f_dim[0], f_dim[1])));
	img_1000_rot.copyTo(win_mat_rot(cv::Rect( f_dim[0]*2, 0, f_dim[0], f_dim[1])));

	cv::namedWindow("Images_rotated");
	cv::imshow("Images_rotated", win_mat_rot);
	cv::waitKey(3000);













	delete [] f_data; f_data = NULL;

	delete [] f_data_double; f_data_double = NULL;
	delete [] rotated_img_0; rotated_img_0 = NULL;
	delete [] imgs_rotated; imgs_rotated = NULL;

}
