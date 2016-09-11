/*
 * PPCA_EM.h
 *
 *  Created on: Sep 11, 2016
 *      Author: kotarokelley
 */

#ifndef PPCA_H_
#define PPCA_H_

class PPCA {
	/** Probabilistic principal component analyzer virtual class. The class and methods are based off of probabilistic PCA model
	 * from "Mixtures of Probabilistic Principal Component Analysers," Tipping and Bishop 1999.
	 *
     *	Class Members:
     *		n_components: integer
     *			 Number of components to keep. If n_components is not set, all
     *			components are kept.
     *		components_: pointer to array of floats
     *			Principal axes in feature space, representing the directions of maximum
     *			variance in the data.
     *		explained_variance_ratio_: pointer to array of floats
     *			Percentage of variance explained by each of the selected components.
     *			Basically eigenvalues normalized by n_components.
     *			Dynamically allocated after calling function fit. Released by deconstructor.
     *		mean_: pointer to array of floats
     *			Per feature empirical mean, estimated from the training set.
     *		noise_variance_: pointer to float
     *			Estimated noise variance.
     *		data: pointer to array of floats.
     *			Dynamically allocated at class instantiation. Released by deconstructor.
     *			     *
     *	Class Methods:
     *
     *		fit:
     *
     *		fit_transform:
     *
     *		get_covariance:
     *
     *		get_params:
     *			Get parameters for this estimator.
     *		get_covariance:
     *			Get covariance of data.
     *		get_covariance_model:
     *			Get covariance of model.
     *
     *
     *
     *
	 */

	public:
		virtual ~ PPCA() {		// clean up dynamically allocated data array
			delete [] data; data = NULL;
			delete [] components_; components_ = NULL;
		}

		virtual void initialize(void) = 0;

		float * components_;

	protected:
		float * data;

};

class PPCA_mixture_EM: public PPCA {

};








#endif /* PPCA_H_ */
