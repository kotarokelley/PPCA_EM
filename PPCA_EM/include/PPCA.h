/*
 * PPCA.h
 *
 *  Created on: Sep 11, 2016
 *      Author: kotarokelley
 */

#ifndef PPCA_H_
#define PPCA_H_

#include "kMeans.h"
#include "Exceptions.h"
#include "Parser.h"
#include "ImageOperations.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <sstream>
#include <armadillo>
#include <limits>

#define SMALL_NUMBER 10e-100
#define LARGE_NUMBER std::numeric_limits<double>::max()

//#include "NumericalRecipes.h"

//using namespace std;
using namespace arma;

class PPCA {
	/** Probabilistic principal component analyzer virtual class. The class and methods are based off of probabilistic PCA model
	 * from "Mixtures of Probabilistic Principal Component Analysers," Tipping and Bishop 1999.
	 *
     *	Class Members:
     *		n_obs: int
     *			Number of observations.
     *		n_var: int
     *			Number of variables for each observation.
     *		n_components: int
     *			Number of components to keep. If n_components is not set, all
     *			components are kept.
     *		mean_: mat
     *			Per feature empirical mean, estimated from the training set.
     *		noise_var: vector of doubles
     *			Estimated noise variance.
     *		data: mat
     *			Input data. Assumes the following format: columns of samples, rows of data points.
     *	Class Methods:
     *		normalize_data:
     *			Subtract mean and divide by std dev for each obs for each variable.
     *		get_n_components
     *			Get how many components we are keeping.
     *		n_models: int
     *			Number of PCA models.
     *		get_mean:
     *			Get mean parameters.
     *		get_noise_variance:
     *			Get estimates for the noise variance of the data.
     *		get_noise_variance_model:
     *			Get estimates for the noise variances of the model.
     *
	 */

	public:
		virtual ~ PPCA() {		// clean up dynamically allocated data array

		}
		/**---Constructors---**/

		PPCA(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models);

		/**---Class Methods---***/

		void normalize_data(void);

		int get_n_obs(void);
			/** Returns the number of observations in data.
			 */
		int get_n_var(void);
			/** Returns the number of variables in each observation.
			 */

	    int get_n_components(void);
			/** Return the number of components that are kept.
			 */
	    int get_n_models(void);
	    	/** Return the number of models in this PCA.
	    	 */
	    mat get_mean(void);
    		/** Return a copy of mean
    		 */
	    std::vector<double> get_noise_variance(void);
    		/** Return a copy of noise_variance
    		 */
	    std::vector<int> get_data_dim(void);
	    	/** Return the dimensions of the indivisual images. xpix, ypix
	    	 */

		/**---Class Members---**/

	    int n_obs;

	    int n_var;

		int n_components;

		int n_models;

		mat mean;

		std::vector<double> noise_var;

		mat data;

		std::vector<int> data_dim;

};

class PPCA_Mixture_EM: public PPCA {
	/** Mixture model probabilistic principal component analysis solved by expectation maximization as derived in
	 *  "Mixtures of Probabilistic Principal Component Analysers," Tipping and Bishop 1999.
	 *
     *	Class Methods:
     *		initialize_random:
     *			Find initial values for each model by assigning by randomly picking data points.
     *		initialize_kmeans:
     *			Find initial values for each model by assigning by kmeans clustering of data points.
     *		optimize:
     *			Execute EM algorithm.
     *		calc_log_Ptn_i:
     *			Helper function to calculate the marginal distribution of a data point for all models.
     *
     *	Class Members:
     *		n_iter: int
     *			Number of iterations so far.
     *		max_iter: int
     *			Maximum number of iterations.
     *		mixfrac: vector of floats
     *			Mixing fraction of each PCA model.
     *		Rni: mat
     *			Posterior responsibility of mixture i for generating data point tn.
     *		W_mat_vector: vector of Mat
     *			Local covariance matrix for each model.
     *		Si_mat_vector: vector of Mat; each matrix represented as 1-D array
     *			Local responsibility weighted covariance matrices for each model.
     *		log_likelihood: double
     *			Value of likelihood function at current iteration.
     *		has_init: int
     *			Flag set to 1 if initialization of all parameters before executing EM. Otherwise 0.
     *
     **/
	public:
		/**---Destructor---**/
		~ PPCA_Mixture_EM() {		// clean up dynamically allocated data array
			//delete [] mixfrac; mixfrac = NULL;

		}

		/**---Constructors---**/

		PPCA_Mixture_EM(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models);
			/**	Multiple models, keep only largest components.
			 */

		PPCA_Mixture_EM(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, char* prev_run_filename);
			/** Use this constructor if starting from a previous job.
			 */

		/**---Class Methods---***/

		int write_to_file_params(int n_iter);  // TODO: for the following write functions, need a standard..
	    	/** Write to file num particles, xpix, ypix, num models, num components,
	    	 * Rni for each particle, Mix frac for each model, mean for each model, Wmat for each model as flat array, Smat for each
	    	 * model as flat array, Minv mat for each model as flat array.
	    	 */
	    void parse_prev_params(char* prev_run_filename);
	    /** Parse saved parameters from previous run.
	     * 	This file should is output by write_to_file_params. The file format is hardcoded.
	     */
	    int write_to_file_mean(int n_iter);
    	/** Write to file mean images as an mrc file.
    	 * 	This is an optional output of the mean for visualization.
    	 */
		void initialize_random(void);
			/** Find initial values for posterior responsibility of mixture i,
			 * 		Rni=p(i|tn) for each model by random clustering of data.
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void initialize_kmeans(int f_max_iter);
			/** Find initial values for posterior responsibility of mixture i,
			 * 		Rni=p(i|tn) for each model by kmeans clustering of data.
			 * 	Arguments:
			 * 		void
			 * 	Returns:
			 * 		void
			 */
		void initialize_helper(void);
			/**	Excecuted by the various initialize functions above to perform the 0th iteration.
			 *
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void update_Rni_all(void);
			/**	Helper function to calculate the posterior responsibility of model i.			 *
			 *	Arguments:
			 *		None.
			 *		i: int
			 *			index for model i.
			 */
		double calc_log_Ptn_i(int f_n, int f_i, mat f_Cinv, double f_log_det_C);
			/**	Helper function to calculate the marginal distribution of a data point for all models.
			 *	Arguments:
			 *		n: int
			 *			index for data point
			 *		i: int
			 *			index for model
			 *		f_Cinv: mat
			 *			pass in pre-calculated Cinv matrix
			 *		f_log_det_C: double
			 *			log determinant of C
			 *	Returns:
			 *		double
			 */
		double calc_log_Ptn_i(mat &f_samples, int f_n, int f_i, mat f_Cinv, double f_log_det_C);
			/**	Helper function to calculate the marginal distribution of a data point for all models.
			 *	Arguments:
			 *		f_samples: mat
			 *			A subset of the data.
			 *		n: int
			 *			index for data point
			 *		i: int
			 *			index for model
			 *		f_Cinv: mat
			 *			pass in pre-calculated Cinv matrix
			 *		f_log_det_C: double
			 *			log determinant of C
			 *	Returns:
			 *		double
			 */
		double calc_log_Ptn_i(mat &f_samples, mat &f_mean, int f_n, int f_i, mat f_Cinv, double f_log_det_C);
			/**	Helper function to calculate the marginal distribution of a data point for all models.
			 *	Arguments:
			 *		f_samples: mat
			 *			A subset of the data.
			 *		f_mean: mat
			 *			Use a different mean than the one store internally for class.
			 *		n: int
			 *			index for data point
			 *		i: int
			 *			index for model
			 *		f_Cinv: mat
			 *			pass in pre-calculated Cinv matrix
			 *		f_log_det_C: double
			 *			log determinant of C
			 *	Returns:
			 *		double
			 */

		void optimize(int f_max_iter,int write_freq_mean, int write_freq_params);
			/** Execute EM algorithm.
			 *	Arguments:
			 *		f_max_iter: int
			 *			Maximum number of iterations to perform.
			 *		write_freq_mean: int
			 *			After how many iterations should we write mean to file. The initial means are always written to file by default.
			 *		write_freq_params: int
			 *			After how many iterations should we write mean to file. The initial params are always written to file by default.
			 *	Returns:
			 *		void
			 */

		/**---Class Members--**/

		std::vector<double> mixfrac;

		int n_iter;

		mat Rni;

		std::vector<mat> W_mat_vector;

		std::vector<mat> S_mat_vector;

		std::vector<mat> Minv_mat_vector;

		double log_likelihood;

		int has_init;


};

class PPCA_Mixture_SAG: public PPCA_Mixture_EM {
	/** Mixture model probabilistic principal component analysis as derived in "Mixtures of Probabilistic Principal Component Analysers," Tipping and Bishop 1999.
	 *  solved by stochastic average gradient method as reported in "A stochastic Gradient Method with an Exponential Convergence Rate for
	 *  Strongly-Convex Optimization with Finite Training Sets," Roux, Schmidt, and Bach, 2012.s
	 *
     *	Class Methods:
     *		initialize_random_SAG:
     *			Find initial values for each model by assigning Rni by randomly picking data points.
     *		initialize_kmeans_SAG:
     *			Find initial values for each model by assigning Rni by kmeans clustering of data points.
     *		get_mini_batch_SAG:
     *			Selects the particles to estimate the gradient at the current iteration.
     *		estimate_Lip_SAG:
     *			Estimate the Lipschitz constant for current sub-level.
     *		calc_base_rate_SAG:
     *			Calculate the base learning rate.
     *		calc_grad_SAG:
     *			Calculates the gradient of the total likelihood.
     *		optimize_SAG:
     *			Execute EM algorithm.
     *	Class Members:
     *		mixfrac_softmax_coeff: std::vector<double>
     *			Coefficients for the softmax coefficients assigned to each mixfrac.
     *		n_mini_batch:	int
     *			size of mini batch.
     *		lip_const: double
     *			estimated Lipschitz constant.
     *		mean_prev: mat
     *			Means from last iteration.
     *		noise_var_prev: std::vector<double>
     *			Noise variance for each model from previous iteration.
     *		mixfrac_prev: std::vector<double>
     *			Mixing fraction of each PCA model from previous iteration.
     *		W_mat_vector_prev: std::vector<mat>
     *			Local covariance matrix for each model from previous iteration.
     *		grad_mean: mat
     *			Running total of gradients of likelihood with respect to the means.
     *		grad_noise_var: std::vector<double>
     *			Running total of gradients of likelihood with respect to the noise variances.
     *		grad_mixfrac: std::vector<double>
     *			Running total of gradients of likelihood with respect to the mix fractions.
     *		grad_W_mat_vector: std::vector<mat>
     *			Running total of gradients of likelihood with respect to factor loadings.
     *		base_rate: double
     *			baseline step size.
     *
     **/
public:
	/**---Destructor---**/
	~ PPCA_Mixture_SAG() {		// clean up dynamically allocated data array
		//delete [] mixfrac; mixfrac = NULL;

	}

	/**---Constructors---**/

	PPCA_Mixture_SAG(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, int f_n_mini_batch);
		/**	Multiple models, keep only largest components.
		 */

	PPCA_Mixture_SAG(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, char* prev_run_filename, int f_n_mini_batch);
		/** Use this constructor if starting from a previous job.
		 */

	/**---Class Methods---***/

	void initialize_random_SAG(void);
		/** Find initial values for posterior responsibility of mixture i,
		 * 		Rni=p(i|tn) for each model by random clustering of data.
		 *	Arguments:
		 *		void
		 *	Returns:
		 *		void
		 */
	void initialize_kmeans_SAG(int f_max_iter);
		/** Find initial values for posterior responsibility of mixture i,
		 * 		Rni=p(i|tn) for each model by kmeans clustering of data.
		 * 	Arguments:
		 * 		void
		 * 	Returns:
		 * 		void
		 */
	std::vector<double> mixfrac_to_softmax(std::vector<double> f_mix_frac_coeff);
	/** Convert from mix_frac vales to their corresponding soft_max coefficients.
	 * 	Arguments:
	 * 		f_mix_frac_coeff: std::vector<double>
	 * 	Returns:
	 * 		std::vector<double>
		 * 			Calculated softmax values.
	 */
	std::vector<double> softmax_to_mixfrac(std::vector<double> f_soft_max_coeff);
		/** Convert from soft_max coefficients to their corresponding mix_frac vales.
		 * 	Arguments:
		 * 		f_soft_max_coeff: std::vector<double>
		 * 	Returns:
		 * 		std::vector<double>
		 * 			Calculated mix_frac values.
		 */
	mat get_mini_batch_SAG(int f_n_mini_batch, double max_trans);
		/**	Randomly selects indices from the unaugmented data and applies a random transformation to each particle.
		 * 	Arguments:
		 * 		f_n_mini_batch:	int
		 * 			size of mini batch
		 * 		max_trans: double
		 * 			maximum translation allowed in units of pixels.
		 * 	Returns:
		 * 		mat:
		 * 			mini batch data.
		 */
	double calc_base_rate_SAG(void);
		/**	Calculate the base learning rate. Return max(1/16, 2^(1-(n_iter/150))): this heruristic is used in Brubaker, Punjani, Fleet, 2015.
		 *	Arguments:
		 *		void
		 *	Returns:
		 *		double:
		 *		base learning rate.
		 */
	std::vector<double> calc_Cinv_log_det_C(std::vector<mat> &f_W_mat_vector, std::vector<double> f_noise_var, std::vector<mat> &f_out_mat_vector);
		/**	Calculate Cinv matrices and corresponding log_det_C values.
		 * 	Arguments:
		 * 		f_W_mat_vector: reference to std::vector<mat>
		 * 			W matrices to use.
		 * 		f_noise_var: std::vector<double>
		 * 			noise_var to use.
		 * 		f_out_mat_vector: reference to std::vector<mat>
		 * 			Output Cinv matrices.
		 *	Returns:
		 *		std::vector<double>
		 *			Calculated values for log_det_C.
		 */
	void calc_log_Ptn_i_mat(std::vector<mat> &f_Cinv_vector, std::vector<double> f_log_det_C_vector, mat &f_mean, mat &f_samples, int f_n_samples, mat &f_out_mat);
		/** Calculate the matrix holding values of log_Ptn_i for a given data set.
		 * 	Arguments:
		 * 		f_Cinv_vector: reference to std::vector<mat>
		 * 			Pre-calculated Cinv matrices
		 * 		f_log_det_C_vector: std::vector<double>
		 * 			Pre-calculated values for log_det_C for each model.
		 * 		f_mean: reference to mat
		 * 			Means used for calculation. n_rows have to match the f_n_samples n_rows and n_cols have to match n_models and will not be checked.
		 * 		f_samples: reference to mat
		 * 			Selected data points to use for calculation.
		 * 		f_n_samples: int
		 * 			Number of samples in f_samples
		 * 		f_out_mat: reference to mat
		 * 			Output matrix. Dimensions should be consisted with data and number of models. Will not be checked.
		 * 	Returns:
		 * 		void
		 *
		 */
	double estimate_log_likelihood(mat &f_log_Ptn_i_mat, std::vector<double> f_mixfrac, mat &f_samples, int f_n_samples);
		/**	Estimate the log likelihood using a subset of the data.
		 * 	Arguments:
		 * 		f_log_Ptn_i_mat: reference to mat
		 * 			Pre-calculated values for Ptn_i.
		 * 		f_mixfrac: std::vector<double>
		 * 			Mix fraction values to use for each model.
		 * 		f_samples: mat
		 * 			A sample of the data set.
		 * 		f_n_samples: int
		 * 			Number of samples .
		 * 	Returns:
		 * 		double:
		 * 			Estimated log likelihood.
		 */

	double estimate_Lip_SAG(int iter, double f_k_1, double f_k, double f_grad_k);
		/** Estimate the Lipschitz constant for current iteration/sub-level by iterative line search.
		 * 		The line search is performed when mod(n_iter,iter) == 0;
		 * 		Else, return n_iter/(2^(1/150)): this heruristic is used in Brubaker, Punjani, Fleet, 2015.
		 * 	Arguments:
		 * 		iter: int
		 * 			Frequency for performing line search.
		 * 		f_samples: reference to mat
		 * 			Sample data to use for use in line search.
		 * 		f_n_samples: int
		 * 			Number of samples.
		 * 	Returns:
		 * 		double:
		 * 			Estimated Lipschitz constant.
		 */

	std::vector<double> calc_grad_SAG(std::vector<mat> &f_Cinv_vector, mat &f_log_Ptn_i_mat,
			std::vector<mat> &f_W_mat_vector, mat &f_mean, std::vector<double> f_mix_frac,
			std::vector<double> f_mix_frac_softmax_coeff, std::vector<double> f_noise_var, mat &f_samples, int f_n_samples );
		/**	Calculate the gradients for parameters by evaluating analytical derivatives.
		 * 	Arguments:
		 * 		f_Cinv_vector: reference to std::vector<mat>
		 * 			Pre-calculated Cinv matrices.
		 * 		f_log_Ptn_i_mat: reference to mat
		 * 			Pre-calculated log_Ptn_i values.
		 * 		f_W_mat_vector:	reference to std::vector<mat>
		 * 			W matrices to use.
		 * 		f_mean:	reference to mat
		 * 			Means to use.
		 * 		f_mix_frac:	std::vector<double>
		 * 			Mix fractions to use.
		 * 		f_mix_frac_softmax_coeff: std::vector<double>
		 * 			Soft-max coefficients used for calculating mix_frac.
		 * 		f_samples: reference to mat
		 * 			Sample data to use for use in line search.
		 * 		f_n_samples: int
		 * 			Number of samples.
		 * 	Returns:
		 * 		std::vector<double>
		 * 			A concatenated vector of all gradient values in the following order: mean(column major), noise_var, mix_frac_softmax, W mat(column major).
		 */
	std::vector<double> calc_grad_finite_dif_SAG(mat &f_log_Ptn_i_mat, std::vector<mat> f_W_mat_vector, mat f_mean, std::vector<double> f_mix_frac,
			std::vector<double> f_mix_frac_softmax_coeff, std::vector<double> f_noise_var, mat &f_samples, int f_n_samples );
	/**	Calculate the gradients for parameters by finite difference method.
	 *
	 */


	/**---Class Members--**/

	std::vector<double> mixfrac_softmax_coeff;

	int n_mini_batch;

	double lip_const;

	double base_rate;

	mat grad_mean;

	std::vector<double> grad_noise_var;

	std::vector<double> grad_mixfrac_softmax_coeff;

	std::vector<mat> grad_W_mat_vector;

	mat grad_mean_prev;

	std::vector<double> grad_noise_var_prev;

	std::vector<double> grad_mixfrac_softmax_coeff_prev;

	std::vector<mat> grad_W_mat_vector_prev;

};


#endif /* PPCA_H_ */
