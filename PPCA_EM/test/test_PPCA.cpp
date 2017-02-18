/*
 * test_PPCA.cpp
 *
 *  Created on: Nov 7, 2016
 *      Author: kotarokelley
 */


#include "Parser.h"
#include "PPCA.h"
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

	std::cout << "Casting data array to double array.\n";

	double * f_data_double = new double[f_dim[0]*f_dim[1]*f_dim[2]];		// convert data to d.

	std::memcpy(f_data_double, f_data, sizeof(double)*f_dim[0]*f_dim[1]*f_dim[2]);

	std::cout << "Constructing a PPCA_Mixture_EM object, initialize with data from: " << filename << "\n";

	std::cout << "Explicitly declaring 2 components and 2 models through constructor.\n";

	arma::mat f_data_mat = arma::mat(f_data_double,f_dim[0]*f_dim[1],f_dim[2]);		// format data into an mat object
	//int f_dim_2 [3];
	//f_dim_2[0] = f_dim[0]; f_dim_2[1] = f_dim[1]; f_dim_2[2] = f_dim[2];

	PPCA_Mixture_EM pca = PPCA_Mixture_EM( f_data_mat, f_dim, 2, 2);

	std::cout << "Checking that data did not get altered by passing to PPCA_Mixture_EM constructor.\n";

	int n_rows = pca.data.n_rows;

	int n_cols = pca.data.n_cols;

	bool equal = true;

	for (int j=0; j<n_cols; j++) {

		for (int i=0; i<n_rows; i++) {

			if (pca.data(i,j) != f_data_double[j*n_rows + i]){

				equal = false;

				break;
			}
		}
	}

	if (equal)
		std::cout << "Data has not been augmented.\n";

	else
		std::cout << "Data has been augmented.\n\n";

	std::clock_t start, stop;

	double elapsed;



	std::cout << "Testing function initialize_uniform.\n";

	start = std::clock();

	pca.initialize_uniform();

	stop = std:: clock();

	elapsed = (stop-start)/(double)CLOCKS_PER_SEC;     //(double(stop-start))/CLOCKS_PER_SEC;

	std::cout << "Function initialize_uniform took: " << elapsed << " s\n\n";


	std::cout << "Displaying the values of Rni for the first image.\n";

	for (int i=0; i<pca.n_models; i++)
		std::cout << pca.Rni(i,0) << "\n";

	std::cout << "\nMake sure that they sum to 1.\n";

	double total = 0.0;

	for (int i=0; i<pca.n_models; i++){

		total += pca.Rni(i,0);
	}

	if (std::abs(total-1.0) < 0.00001)
		std::cout << "Probabilities Rni for first images summed to 1.\nPass!\n\n";

	else
		std::cout << "Probabilities Rni for first images did not sum to 1.\n" << "Summed to: " << total << "\nNoPass\n\n";

	std::cout << "Testing function write_to_file\n";

	pca.write_to_file_params("test_PPCA_parameters.bmp");

	std::cout << "\nParameters have been saved to file: test_PPCA_parameters.bmp\n";

	std::cout << "\nLets check to see if the data didn't get corrupted during writing\n";

	std::cout << "\nConstructing a new pca object taking in parameters from test_PPCA_parameters.bmp";

	PPCA_Mixture_EM pca2 = PPCA_Mixture_EM( f_data_mat, f_dim, 2, 2, "test_PPCA_parameters.bmp");


	if (pca.n_obs != pca2.n_obs)
		std::cout << "\nb_obs not consistent.\n";

	else
		std::cout << "\nb_obs is consistent.\n";


	if (pca.n_var != pca2.n_var)
		std::cout << "\nn_var not consistent.\n";

	else
		std::cout << "\nn_var is consistent.\n";

	if (pca.n_components != pca2.n_components)
		std::cout << "\nn_components not consistent.\n";

	else
		std::cout << "\nn_components is consistent.\n";

	if (pca.n_models != pca2.n_models)
		std::cout << "\nn_models not consistent.\n";

	else
		std::cout << "\nn_models is consistent.\n";

	equal = 1;
	for(int i=0; i<pca.n_obs*pca.n_models; i++){
		if(pca.Rni(i) != pca2.Rni(i))
			equal = 0;
	}
	if (!equal)
			std::cout << "\nRni not consistent\n";
		else
			std::cout << "\nRni consistent\n";

	equal = 1;
	for (int i=0; i<pca.n_models; i++){

		if (pca.mixfrac[i] != pca2.mixfrac[i])
			equal = 0;
	}

	if (!equal)
		std::cout << "\nmixfrac not consistent\n";
	else
		std::cout << "\nmixfrac consistent\n";

	equal = 1;
	for (int i=0; i<pca.n_models*pca.n_var; i++){
		if (pca.mean[i] != pca2.mean[i])
			equal = 0;
	}

	if (!equal)
		std::cout << "\nmean not consistent\n";
	else
		std::cout << "\nmean consistent\n";

	equal = 1;
	for (int i=0; i<pca.n_models; i++){
		for (int j=0; j<pca.n_var*pca.n_components; j++){
			if (pca.W_mat_vector[i](j) != pca2.W_mat_vector[i](j))
				equal = 0;
		}
	}

	if (!equal)
		std::cout << "\nWmat not consistent\n";
	else
		std::cout << "\nWmat consistent\n";

	equal = 1;
	for (int i=0; i<pca.n_models; i++){
		for (int j=0; j<pca.n_var*pca.n_var; j++){
			if (pca.S_mat_vector[i](j) != pca2.S_mat_vector[i](j))
				equal = 0;
		}
	}

	if (!equal)
		std::cout << "\nSmat not consistent\n";
	else
		std::cout << "\nSmat consistent\n";

	equal = 1;
	for (int i=0; i<pca.n_models; i++){
		for (int j=0; j<pca.n_components*pca.n_components; j++){
			if (pca.Minv_mat_vector[i](j) != pca2.Minv_mat_vector[i](j))
				equal = 0;
		}
	}

	if (!equal)
		std::cout << "\nMinv_mat not consistent\n";
	else
		std::cout << "\nMinv_mat consistent\n";


	std::cout << "Testing function initialize_random.\n";
	start = std::clock();

	pca.initialize_random();
	stop = std::clock();
	elapsed = (double(stop-start))/CLOCKS_PER_SEC;
	std::cout << "Function initialize_random took: " << elapsed << " s\n\n";

	std::cout << "Displaying the values of Rni for the first image.\n";
	for (int i=0; i<pca.n_models; i++){
		std::cout << pca.Rni(0,i) << "\n";
	}
	std::cout << "\nMake sure that they sum to 1.\n";

	total = 0.0;

	for (int i=0; i<pca.n_models; i++){
		total += pca.Rni(0,i);
	}
	if (std::abs(total-1.0) < 0.00001)
		std::cout << "\nProbabilities Rni for first images summed to 1.\nPass!\n";
	else
		std::cout << "\nProbabilities Rni for first images did not sum to 1.\n" << "Summed to: " << total << "\nNoPass\n";


	std::cout << "Performing one round of optimizations.\n";
	start = std::clock();

	pca.optimize(2);
	stop = std::clock();
	elapsed = (double(stop-start))/CLOCKS_PER_SEC;
	std::cout << "Function optimize took: " << elapsed << " s\n\n";
	std::cout << "Writing to file Rni\n";







	// cleanup
	//delete testParser; testParser = NULL;
	delete f_data; f_data = NULL;
	delete f_header; f_header = NULL;
	delete [] f_data_double; f_data_double = NULL;
	//delete pca; pca = NULL;

}




