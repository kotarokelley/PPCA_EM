/*
 * test_kMeans.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: kotarokelley
 */




#include "kMeans.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <random>
#include <fstream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui.hpp"

using namespace arma;

fmat init_board_gaus(int N, int k){
	/** Returns fmat consisting of columns of data points coming from k centers.
	 * 	Each data point has dimension of 2 (2 rows).
	 *
	 */
	int n = N/k;
	fmat out_data(2,N);
	std::srand(std::time(NULL));										// generate seed.
	std::default_random_engine generator;

	for (int i=0; i<k; i++){
		float ux = -1 + static_cast<float> (std::rand()/(static_cast<float>(RAND_MAX/(1+1)))); // generate mean between -1 and 1
		float uy = -1 + static_cast<float> (std::rand()/(static_cast<float>(RAND_MAX/(1+1))));
		float sig = 0.05 + static_cast<float> (std::rand()/(static_cast<float>(RAND_MAX/(0.5-0.05))));

		std::normal_distribution<float> norm_x(ux, sig);
		std::normal_distribution<float> norm_y(uy, sig);

		int counter = 0;
		while(counter<n){
			float dat_x = norm_x(generator);
			float dat_y = norm_y(generator);
			if (dat_x<1 && dat_y<1){
				out_data(0,n*i + counter) = dat_x;
				out_data(1,n*i + counter) = dat_y;
				counter++;

			}
		}


	}

	return out_data;
}


//using namespace arma;

int main(void){

	std::cout << "Testing class kMeans.\n\n";
	std::cout << "Generating 200 data points from 3 gaussian distributions.\n";

	fmat data = init_board_gaus(200,3);					// Generate data.

	float * data_x = new float[data.n_cols];			// Convert data to simple array for feeding into opencv
	float * data_y = new float[data.n_cols];
	for (int i=0; i<data.n_cols; i++){
		data_x[i] = data(0,i);
		data_y[i] = data(1,i);
	}
	std::ofstream file;									// Open file to write data for analysis with python script.
	file.open("file.txt");
	file<< "x_data:\n";
	for (int i=0; i<data.n_cols; i++){
		file<< data_x[i];
		file<< " ";
	}
	file<< "\ny_data:\n";
	for (int i=0; i<data.n_cols; i++){
		file<< data_y[i];
		file<< " ";
	}
	file<< "\n";

	kMeans classifier(3,20,data);
	classifier.initCentersRandom();
	std::vector<int> solution;
	solution = classifier.findSolution();

	file<< "kMeans_solution\n";
	for (int i=0; i<data.n_cols; i++){
		file<< solution[i];
		file<< " ";
	}
	file<< "\n";





	file.close();
	delete data_x; data_x = NULL;
	delete data_y; data_y = NULL;
}
