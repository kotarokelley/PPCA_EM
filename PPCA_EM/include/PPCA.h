/*
 * PPCA.h
 *
 *  Created on: Sep 11, 2016
 *      Author: kotarokelley
 */

#ifndef PPCA_H_
#define PPCA_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <armadillo>
#include <vector>
//#include "NumericalRecipes.h"

//using namespace std;
using namespace arma;

class PPCA {
	/** Probabilistic principal component analyzer virtual class. The class and methods are based off of probabilistic PCA model
	 * from "Mixtures of Probabilistic Principal Component Analysers," Tipping and Bishop 1999.
	 *
     *	Class Members:
     *		n_components: int
     *			Number of components to keep. If n_components is not set, all
     *			components are kept.
     *		components_: array of fmat
     *			Principal axes in feature space, representing the directions of maximum
     *			variance in the data.
     *		explained_variance_ratio_: array of floats
     *			Percentage of variance explained by each of the selected components.
     *			Basically eigenvalues normalized by n_components.
     *			Dynamically allocated after calling function fit. Released by deconstructor.
     *		mean_: array of frowvec
     *			Per feature empirical mean, estimated from the training set.
     *		noise_variance_: array of frowvec
     *			Estimated noise variance.
     *		noise_variance_model: array of frowvec
     *			Estimated noise variance for model.
     *		data: fmat
     *			Input data
     *	Class Methods:
     *
     *		fit:
     *			Fit the data to the model.
     *		fit_transform:
     *			Fit the data and apply the dimensionality reduction on data. TODO >>> Decide if passed in data is modified.
     *		get_covariance:
     *			Get covariance of data.
     *		get_covariance_model:
     *			Get covariance of model.
     *		get_params:
     *			Get parameters for this estimator.
	 */

	public:
		virtual ~ PPCA() {		// clean up dynamically allocated data array

		}

		/**---Class Methods---***/

		virtual void fit(void) = 0;
			/** Top level function called to perform PCA analysis.
			 *
			 */
	    virtual char * get_params(void) = 0;
			/**	Arguments:
			 * 		void
			 * 	Returns:
			 * 		String representing the parameters for this model.
			 */
	    virtual int get_n_components(void) = 0;
			/** Return a copy of n_componentes.
			 */
	    virtual std::vector<fmat> get_components_(void) = 0;
	    	/** Return a copy of componentes_.
	    	 */
	    virtual std::vector<frowvec> get_explained_variance_ratio(void) = 0;
    		/** Return a copy of explained_variance_ratio.
    		 */
	    virtual std::vector<frowvec> get_mean(void) = 0;
    		/** Return a copy of mean
    		 */
	    virtual std::vector<frowvec> get_noise_variance(void) = 0;
    		/** Return a copy of noise_variance
    		 */
	    virtual std::vector<frowvec> get_noise_variance_moedl(void) = 0;
    		/** Return a copy of noise_variance_model
    		 */

		/** Class Members **/
	protected:
		int n_components;

		std::vector<fmat> components_;

		std::vector<frowvec> explained_variance_ratio;

		std::vector<frowvec> mean;

		std::vector<frowvec> noise_variance;

		std::vector<frowvec> noise_variance_model;

		fmat data;

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
     *		n_models: int
     *			Number of PCA models.
     *		mix_frac_array: array of floats
     *			Mixing fraction of each PCA model.
     *		W_mat_array: array of Mat
     *			Local covariance matrix for each model.
     *		Si_mat_array: array of Mat; each matrix represented as 1-D array
     *			Local responsibility weighted covariance matrices for each model.
     **/
	public:
		/**---Destructor---**/
		~ PPCA_Mixture_EM() {		// clean up dynamically allocated data array
			delete [] mixfrac; mixfrac = NULL;
			delete [] Rni; Rni = NULL;
			delete [] Rn; Rn = NULL;
		}
		/**---Constructors---**/
		PPCA_Mixture_EM(fmat f_data, int* f_dim);
			/**	Single model, keep all components.
			 */
		PPCA_Mixture_EM(fmat f_data, int* f_dim, int f_n_components);
			/**	Single model, keep only largest components.
			 */
		PPCA_Mixture_EM(fmat f_data, int* f_dim, int f_n_components, int f_n_models);
			/**	Multiple models, keep only largest components.
			 */
		/**---Class Methods---***/

		void fit(void);

	    char * get_params(void);
			/**	Arguments:
			 * 		void
			 * 	Returns:
			 * 		String representing the parameters for this model.
			 **/
		int get_n_components(void) = 0;
			/** Return a copy of n_componentes.
			 */
		std::vector<fmat> get_components_(void) = 0;
			/** Return a copy of componentes_.
			 */
		std::vector<frowvec> get_explained_variance_ratio(void) = 0;
			/** Return a copy of explained_variance_ratio.
			 */
		std::vector<frowvec> get_mean(void) = 0;
			/** Return a copy of mean
			 */
		std::vector<frowvec> get_noise_variance(void) = 0;
			/** Return a copy of noise_variance
			 */
		std::vector<frowvec> get_noise_variance_moedl(void) = 0;
			/** Return a copy of noise_variance_model
			 */

		void initialize_random(void);
			/** Find initial values for mixing fraction pi and mean mu for each model by initial clustering of data.
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
		void optimize(float threshold);
			/** Execute EM algorithm.
			 *	Arguments:
			 *		void
			 *	Returns:
			 *		void
			 */
			/**---Class Members--**/

		int n_components;

		std::vector<fmat> components_;

		std::vector<frowvec> explained_variance_ratio;

		std::vector<frowvec> mean;

		std::vector<frowvec> noise_variance;

		std::vector<frowvec> noise_variance_model;


		int n_models;

		float * mixfrac;

		float * Rni;

		float * Rn;

		std::vector<fmat> W_mat_vector;

		std::vector<fmat> Si_mat_vector;

};


#endif /* PPCA_H_ */
