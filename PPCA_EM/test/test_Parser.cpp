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

	mrcParser  testParser =  mrcParser(filename);  // input binary data

	char * f_data = testParser.getData(false);			// extract data, verbose off

	char * f_header = testParser.getHdr();

	std::vector<int> f_dim = testParser.getDim();					// get dimension of image stack

	int f_pixelType = testParser.getPixelType();		// get pixel type

	testParser.hdrInfo();

	std::cout << "\nCasting data array to float array.\n\n";
	float * f_data_float = new float[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.
	std::memcpy(f_data_float, f_data, sizeof(float)*f_dim[0]*f_dim[1]*f_dim[2]);

	std::cout << "\nCastering data array from float to double and then back to double.\n";
	double * f_data_double = new double[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.
	std::copy(f_data_float, f_data_float+f_dim[0]*f_dim[1]*f_dim[2], f_data_double);
	std::copy(f_data_double, f_data_double+f_dim[0]*f_dim[1]*f_dim[2], f_data_float);

	mrcParser writeParser = mrcParser("test_Parser.mrc",0);

	std::cout << "\nWriting header to file: test_Parser.mrc\n";
	writeParser.writeHdr();


	cv::Mat img_0(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,f_data_float);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_0);
	cv::waitKey(2000);

	std::cout << "\nWriting data to file: test_Parser.mrc\n";
	writeParser.writeData(f_data_float, f_dim);

	//delete writeParser; writeParser = NULL;						// Close parser. File is closed by destructor.

	mrcParser readParser = mrcParser("test_Parser.mrc");
	readParser.parseHdr();
	readParser.hdrInfo();

	char * r_data = readParser.getData(false);			// extract data, verbose off


	std::cout << "\nCasting data array to float array.\n\n";
	float * r_data_float = new float[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.
	std::memcpy(r_data_float, r_data, sizeof(float)*f_dim[0]*f_dim[1]*f_dim[2]);

	std::cout << "Displaying first image again after transfer.\n";

	cv::Mat img_1(cv::Size(f_dim[0],f_dim[1]),CV_32FC1,r_data_float);
	cv::namedWindow("image_0",cv::WINDOW_NORMAL);
	cv::imshow("image_0", img_1);
	cv::waitKey(2000);									// display for 5 sec to check that image is not augmented.

	delete f_data; f_data = NULL;

	delete f_data_double; f_data_double = NULL;
	delete r_data; r_data = NULL;
	//delete r_data_double; r_data_double = NULL;

}
