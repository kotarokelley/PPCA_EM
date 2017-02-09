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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <tuple>
#include <math.h>
#include <armadillo>


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
     *		components_: array of mat
     *			Principal axes in feature space, representing the directions of maximum
     *			variance in the data.
     *		explained_variance_ratio_: array of doubles
     *			Percentage of variance explained by each of the selected components.
     *			Basically eigenvalues normalized by n_components.
     *			Dynamically allocated after calling function fit. Released by deconstructor.
     *		mean_: mat
     *			Per feature empirical mean, estimated from the training set.
     *		noise_var: vector of doubles
     *			Estimated noise variance.
     *		noise_var_model: vector of doubles
     *			Estimated noise variance for model.
     *		data: mat
     *			Input data. Assumes the following format: columns of samples, rows of data points.
     *	Class Methods:
     *
     *		get_params:
     *			Get parameters for this estimator.
     *		get_n_components_
     *			Get how many components we are keeping.
     *		n_models: int
     *			Number of PCA models.
     *		get_components:
     *			Get a copy of components_
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

		//virtual std::tuple<int, int, int, int> get_params(void) = 0;
			/**	Arguments:
			 * 		void
			 * 	Returns:
			 * 		Tuple representing the parameters for this model.
			 * 			Rows, Cols, n_components, n_models.
			 */
	    virtual int get_n_components(void) = 0;
			/** Return a copy of n_components.
			 */
	    //virtual std::vector<mat> get_components_(void) = 0;
	    	/** Return a copy of components_.
	    	 */
	    virtual mat get_mean(void) = 0;
    		/** Return a copy of mean
    		 */
	    virtual std::vector<double> get_noise_variance(void) = 0;
    		/** Return a copy of noise_variance
    		 */
	    virtual std::vector<int> get_data_dim(void) = 0;
	    	/** Return the dimensions of the indivisual images. xpix, ypix
	    	 */
	    //virtual std::vector<float> get_noise_variance_model(void) = 0;
    		/** Return a copy of noise_variance_model
    		 */

		/** Class Members **/

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
     *			Find initi{al values for mixing fraction and mean for each model by randomly picking data points.
     *		optimize:
     *			Execute EM algorithm.
     *	Class Members:
     *		n_iter: int
     *			Number of iterations so far.
     *		max_iter: int
     *			Maximum number of iterations.
     *		mixfrac: vector of floats
     *			Mixing fraction of each PCA model.
     *		Rni: fmat
     *			Posterior responsibility of mixture i for generating data point tn.
     *		W_mat_vector: vector of Mat
     *			Local covariance matrix for each model.
     *		Si_mat_vector: vector of Mat; each matrix represented as 1-D array
     *			Local responsibility weighted covariance matrices for each model.
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

		std::tuple<int, int, int, int> get_params(void);
			/**	Arguments:
			 * 		void
			 * 	Returns:
			 * 		Tuple representing the parameters for this model.
			 * 			Rows, Cols, n_components, n_models.
			 */
		int get_n_components(void);
			/** Return a copy of n_components.
			 */
		//std::vector<mat> get_components_(void);
			/** Return a copy of componentes_.
			 */
		//std::vector<rowvec> get_explained_variance_ratio(void);
			/** Return a copy of explained_variance_ratio.
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
		int write_to_file_params(char* filename);  // TODO: for the following write functions, need a standard..
	    	/** Write to file num particles, xpix, ypix, num models, num components,
	    	 * Rni for each particle, Mix frac for each model, mean for each model, Wmat for each model as flat array, Smat for each
	    	 * model as flat array, Minv mat for each model as flat array.
	    	 */
	    void parse_prev_params(char* prev_run_filename);
	    /** Parse saved parameters from previous run.
	     * 	This file should is output by write_to_file_params. The file format is hardcoded.
	     */
	    int write_to_file_mean(char* filename);
    	/** Write to file mean images as an mrc file.
    	 * 	This is an optional output of the mean for visualization.
    	 */

		void initialize_uniform(void);
			/** Find initial values for posterior responsibility of mixture i,
			 * 		Rni=p(i|tn) for each model by assigning uniform value of 1/n_models.
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void initialize_random(void);
			/** Find initial values for posterior responsibility of mixture i,
			 * 		Rni=p(i|tn) for each model by random clustering of data.
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void initialize_kmeans(void);
			/** Find initial values for posterior responsibility of mixture i,
			 * 		Rni=p(i|tn) for each model by kmeans clustering of data.
			 * 	Arguments:
			 * 		void
			 * 	Returns:
			 * 		void
			 */
		void initialize_helper(void);
			/*	Excecuted by the various initialize functions above to perform the 0th iteration.
			 *
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void update_Rni(int n);
			/*	Helper function to calculate the posterior responsibility of model i.			 *
			 *	Arguments:
			 *		n: int
			 *			index for data point
			 */
		double calc_Ptn_i(int n, int i);				//TODO this function is missing a constant that drops out in the Rni equation. Fix?
			/*	Helper function to calculate the marginal distribution of a data point for all models.
			 *	Arguments:
			 *		n: int
			 *			index for data point
			 *		i: int
			 *			index for model
			 */
		void optimize(int f_max_iter);
			/** Execute EM algorithm.
			 *	Arguments:
			 *		f_max_iter: int
			 *			Maximum number of iterations to perform.
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

		int has_init;

};


#endif /* PPCA_H_ */
