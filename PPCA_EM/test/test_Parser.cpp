/*
 * test_Parser.cpp
 *
 *  Created on: Feb 6, 2017
 *      Author: kotarokelley
 */

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

	mrcParser * testParser = new mrcParser(filename);  // input binary data

	char * f_data = testParser->getData(false);			// extract data, verbose off

	char * f_header = testParser->getHdr();

	int * f_dim = testParser->getDim();					// get dimension of image stack

	int f_pixelType = testParser->getPixelType();		// get pixel type

	testParser->hdrInfo();

	std::cout << "\nCasting data array to double array.\n\n";
	double * f_data_double = new double[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.
	std::memcpy(f_data_double, f_data, sizeof(double)*f_dim[0]*f_dim[1]*f_dim[2]);

	mrcParser * writeParser = new mrcParser("test_Parser.mrc",0);

	//writeParser->num[0] = f_dim[0];
	//writeParser->num[1] = f_dim[1];
	//writeParser->num[2] = f_dim[2];
	//writeParser->pixelType = f_pixelType;

	std::cout << "\nWriting header to file: test_Parser.mrc\n";
	writeParser->writeHdr();

	std::cout << "\nCasting f_data array to float array.\n\n";
	float * f_data_float = new float[f_dim[0]*f_dim[1]*f_dim[2]];
	std::copy(f_data_double, f_data_double + f_dim[0]*f_dim[1]*f_dim[2], f_data_float);

	std::cout << "\nWriting data to file: test_Parser.mrc\n";
	writeParser->writeData(f_data_float, f_dim);

	delete writeParser; writeParser = NULL;						// Close parser. File is closed by destructor.

	mrcParser * readParser = new mrcParser("test_Parser.mrc");
	readParser->parseHdr();
	readParser->hdrInfo();

	char * r_data = readParser->getData(false);			// extract data, verbose off


	std::cout << "\nCasting data array to double array.\n\n";
	double * r_data_double = new double[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.
	std::memcpy(r_data_double, f_data, sizeof(double)*f_dim[0]*f_dim[1]*f_dim[2]);

	std::cout << "Displaying first image.\n";

	cv::Mat img_0(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,r_data_double);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_0);
	cv::waitKey(5000);									// display for 5 sec to check that image is not augmented.

	delete f_data; f_data = NULL;
	delete f_data_double; f_data_double = NULL;
	delete r_data; r_data = NULL;
	delete r_data_double; r_data_double = NULL;

}
