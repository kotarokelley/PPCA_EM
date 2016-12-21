/*
 * test_class2D.cpp
 *
 *  Created on: Oct 8, 2016
 *      Author: kotarokelley
 */

#include "Parser.h"
#include "class2D.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui.hpp"


int main(void){

	char * filename = "zip_test.mrc";

	Parser * testParser = new mrcParser(filename);  // input binary data

	char * f_data = testParser->getData(false);			// extract data, verbose off

	int * f_dim = testParser->getDim();					// get dimension of image stack

	int f_pixelType = testParser->getPixelType();		// get pixel type


	std::cout << "Data Set:\n";
	std::cout << filename << "\n";
	std::cout << "pixel type\n";
	std::cout << f_pixelType << "\n";

	std::cout << "dimensions of data set\n";
	std::cout << f_dim[0] << "\n";
	std::cout << f_dim[1] << "\n";
	std::cout << f_dim[2] << "\n";

	std::cout << "Casting data array to float array.\n";
	float * f_data_float = new float[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to float.
	std::memcpy(f_data_float, f_data, sizeof(float)*f_dim[0]*f_dim[1]*f_dim[2]);

	std::cout << "Displaying first image.\n";

	cv::Mat img_0(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,f_data_float);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_0);
	cv::waitKey(1000);									// display for 1 sec

	std::cout << "Constructing a class2D object, initialize with data from: " << filename << "\n";
	class2D classifier = class2D(f_data_float, f_dim, 10);
	std::cout << "Checking that data did not get altered by passing to class2D constructor.\n";
	int n_rows = classifier.data.n_rows;
	int n_cols = classifier.data.n_cols;
	bool equal = true;

	for (int j=0; j<n_cols; j++) {
		for (int i=0; i<n_rows; i++) {

			if (classifier.data(i,j) != f_data_float[j*n_rows + i]){
				equal = false;
				break;
			}
		}
	}

	if (equal)
		std::cout << "Data has not been augmented.\n";
	else
		std::cout << "Data has been augmented.\n";

	std::cout << "Displaying first image again.\n";
	float * f_data_float_mod = new float[f_dim[0]*f_dim[1]];

	for (int i=0; i<f_dim[0]*f_dim[1]; i++) {
		f_data_float_mod[i] = classifier.data(i,0);
	}

	cv::Mat img_1(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,f_data_float);
	cv::namedWindow("image_1",cv::WINDOW_NORMAL);
	cv::imshow("image_1", img_1);
	cv::waitKey(1000);

	std::cout << "Testing method class2D::classify_PPCA_EM.\n";
	classifier.classify_PPCA_EM();












	// cleanup
	delete testParser; testParser = NULL;
	delete f_data; f_data = NULL;
	delete [] f_data_float; f_data_float = NULL;
	delete [] f_data_float_mod; f_data_float_mod = NULL;



	/**
	cv::Mat img_0(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,f_data);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_0);
	cv::waitKey(0);
	// Use std::atof to convert from char to 32 bit float.


	cv::Mat ogImg(cv::Size(this->cols,this->rows),CV_32FC1,this->data);
	cv::Mat reImg;
	cv::resize(ogImg,reImg,cv::Size(this->cols*4,this->rows*4));
	double min;
	double max;
	cv::minMaxLoc(reImg,&min,&max);
	cv::Mat adjMap;
	reImg.convertTo(adjMap,CV_8U,255.0/(max-min),-min*225.0/(max-min));
	//cv::convertScaleAbs(reImg,adjMap,255/max);

	cv::namedWindow("image",cv::WINDOW_NORMAL);
	cv::imshow("image",adjMap);
	cv::waitKey(0);
	***/
}
