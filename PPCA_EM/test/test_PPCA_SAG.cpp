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


int main(void){

	char * filename = "zip_test_plusnoise.mrc";

	mrcParser testParser = mrcParser(filename);  					// input binary data

	char * f_data = testParser.getData(false);						// extract data, verbose off

	char * f_header = testParser.getHdr();

	std::vector<int> f_dim = testParser.getDim();					// get dimension of image stack

	int f_pixelType = testParser.getPixelType();					// get pixel type


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

	std::cout << "Constructing a PPCA_Mixture_SAG object, initialize with data from: " << filename << "\n";

	std::cout << "Declare 2 components and 2 models mini batch size: 10 through constructor.\n";

	arma::mat f_data_mat = arma::mat(f_data_double,f_dim[0]*f_dim[1],f_dim[2]);		// format data into an mat object


	PPCA_Mixture_SAG pca = PPCA_Mixture_SAG( f_data_mat, f_dim, 2, 2, 10);

	std::cout << "Checking that data did not get altered by passing to PPCA_Mixture_SAG constructor.\n";

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
	double total = 0.0;


	std::cout << "Testing function initialize_random.\n";
	start = std::clock();
	pca.initialize_random_SAG();

	std::cout << "\n";
	stop = std::clock();
	elapsed = (double(stop-start))/CLOCKS_PER_SEC;
	std::cout << "Function initialize_random took: " << elapsed << " s\n\n";

	std::cout << "Displaying the values of mixfrac for each model.\n";

	total = 0.0;
	for (int i=0; i<pca.n_models; i++){
		std::cout << "mixfrac[" << i << "]: " << pca.mixfrac[i] << "\n";
		total += pca.mixfrac[i];
	}
	if (std::abs(total-1.0) < 0.00001)
		std::cout << "\nProbabilities of mixfrac summed to 1. \nPass!\n";
	else
			std::cout << "\nProbabilities mixfrac did not sum to 1.\n" << "Summed to: " << total << "\nNoPass\n";


	std::cout << "Displaying the values of mixfrac_softmax_coeff for each model.\n\n";

	for (int i=0; i<pca.n_models; i++)
		std::cout << "mixfrac_softmax_coeff[" << i << "]: " << pca.mixfrac_softmax_coeff[i] << "\n";
	std::cout << "\n";

	std::cout << "Computing gradient vector using function calc_grad_SAG\n\n";

	int n_samples(10);
	int n_var = pca.n_var;
	int n_models = pca.n_models;
	std::vector<mat> C_inv_vector(n_models,mat(n_var,n_var,fill::zeros));
	mat log_Ptn_i_mat(n_samples,n_models);

	mat sample_data(n_var,n_samples);
	std::srand(std::time(NULL));							// generate seed.

	for (int n=0; n<n_samples; n++){
		int idx = std::rand() % 100;
		sample_data.col(n) = pca.data.col(idx);
	}

	std::vector<double> log_det_C = pca.calc_Cinv_log_det_C(pca.W_mat_vector, pca.noise_var, C_inv_vector);
	pca.calc_log_Ptn_i_mat(C_inv_vector, log_det_C, pca.mean, sample_data, n_samples, log_Ptn_i_mat);

	std::vector<double> grad_analytical;
	for (int i=0; i<10; i++)
		grad_analytical = pca.calc_grad_SAG(C_inv_vector,log_Ptn_i_mat, pca.W_mat_vector, pca.mean, pca.mixfrac, pca.mixfrac_softmax_coeff, pca.noise_var, sample_data, n_samples);


	std::cout << "Computing gradient vector using function calc_grad_finite_dif_SAG\n\n";

	std::vector<double> grad_finite_dif;
	for (int i=0; i<10; i++)
		grad_finite_dif = pca.calc_grad_finite_dif_SAG(pca.W_mat_vector, pca.mean, pca.mixfrac, pca.mixfrac_softmax_coeff, pca.noise_var, sample_data, n_samples );

	std::cout << "Comparing gradient vectors calculated analytically and finite differences\n\n";

	for (int k=0; k<grad_analytical.size(); k++){
			std::cout << grad_analytical[k] << "	"<< grad_finite_dif[k] << "\n";
	}
	std::cout << "\n";


	//std::cout << "Performing 100 round of optimizations.\n";
	//start = std::clock();

	//pca.optimize(20,1,2);
	//stop = std::clock();
	//elapsed = (double(stop-start))/CLOCKS_PER_SEC;
	//std::cout << "Function optimize took: " << elapsed << " s\n\n";





	// cleanup
	//delete testParser; testParser = NULL;
	//delete f_data; f_data = NULL;
	//delete f_header; f_header = NULL;



	delete [] f_data_float; f_data_float = NULL;
	delete [] f_data_double; f_data_double = NULL;
	//delete pca; pca = NULL;

}








