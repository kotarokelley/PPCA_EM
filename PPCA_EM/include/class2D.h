/*
 * class2D.h
 *
 *  Created on: Oct 8, 2016
 *      Author: kotarokelley
 */

#ifndef CLASS2D_H_
#define CLASS2D_H_
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <armadillo>
#include <vector>
#include "PPCA.h"



using namespace arma;

class class2D{
	/** Classify 2D particles by PPCA.
	 *
	 *	Class Members:
	 *		data: fmat
	 *			Input data rearranged to be columns of samples, rows of data points.
	 *		data_dim: int array size 2
	 *			The dimensions of data. This information can be accessed directly through data fmat object so maybe delete?
	 *		n_classes: int
	 *			Number of different classes.
	 *		n_components: int
	 *			Number of components to keep in PPCA model. Integer between (1,num variables). Smaller
	 *			number will be faster, but larger number will be more accurate.
	 *			If not designated, it will be determined automatically.
	 *

	 *	Class Methods:
	 *		classify:
	 *			Classify 2D particles.
	 *		tofile:
	 *			Output alignment and classification paramters.
	 **/
	public:


	~ class2D() {		// clean up dynamically allocated data array
		//delete [] data; data = NULL;
	}


	/**---Class Constructor---***/
	class2D(float* f_data, int* f_dim, int f_n_classes);
	/**
	 * f_data:
	 * 		pointer to a 1-D array of floats. Internally, this is converted to an arma::fmat data structure.
	 * f_dim:
	 * 		array describing the dimensions of the data set. x pixels, y pixels, n images.
	 */
	class2D(float* f_data, int* f_dim, int f_n_classes, int f_n_components);

	/**---Class Methods---***/
    void classify_PPCA_EM(void);
    /** Classify 2D particles by PPCA.
     */
    void tofile(char * f_name);
    /**Output alignment and classification parameters.
     */

    /**---Class members---**/
    arma::fmat data;
    int data_dim[2];
    int n_classes;
    int n_components;
};




#endif /* CLASS2D_H_ */
