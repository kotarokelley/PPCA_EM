/*
 * kMeans.h
 *
 *  Created on: Nov 13, 2016
 *      Author: kotarokelley
 */

#ifndef KMEANS_H_
#define KMEANS_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <armadillo>

using namespace arma;

class kMeans{
	/** Implementation of Lloyd's algorithm fokMeans clustering. This class is specifically designed to
	 * 	output initial values of Rni for class PPCA_Mixture_EM.
	 *
	 *	Class Members:
	 *		numDat: int
	 *			Number of data points.
	 *		numCenters: int
	 *			Number of centers.
	 *		numIter: int
	 *			Maximum number of iterations before terminations if convergence under delta has not been reached.
	 *		data: fmat
	 *			Input data. Assumes the following format: columns of samples, rows of data points.
	 *		centers: std::vector<float>
	 *			Hold values for the centers.
	 *		newgroup: std::vector<int>
	 *			Indices of the center assigned to each data point from current iteration.
	 *		oldgroup: std::vector<int>
	 *			Indices of the center assigned to each data point from previous iteration.
	 *		Rni_init: fmat
	 *			Initial values for PPCA_Mixture_EM.
	 *	Class Methods:
	 *		initCentRandom:
	 *			Initialize centers by randomly choosing numCenters from the data set.
	 *		findSolutions:
	 *			Run Lloyd's algorithm until convergence under delta or numInter iterations.
	 */


	public:

		/**---Class Constructors---**/
		kMeans(int f_numCenters, int f_numIter, mat f_data);

		/**---Class Methods---***/
		void initCentersRandom(void);
		/** Initialize centers by randomly choosing numCenters from the data set.
		 *	Arguments:
		 *		void
		 *	Returns:
		 *		void
		 */
		std::vector<int> findSolution(void);
		/** Initialize centers by randomly choosing numCenters from the data set.
		 *	Arguments:
		 *		void
		 *	Returns:
		 *		newgroup: std::vector<int>
		 *			Indices of the center assigned to each data point.
		 */
		mat getCenters(void);
		/** Get the means.
		 *
		 */
		int findClosest(int in_indx);
		/**	Find closest center for input inVec.
		 * 	Arguments:
		 * 		in_indx: int
		 * 			Column of this->data to analyze.
		 * 	Returns:
		 * 		Index of closest center.
		 */
		void recalc_Centers(void);
		/**	Find closest center for input inVec.
		 * 	Arguments:
		 *
		 * 	Returns:
		 */
		bool hasConverged(void);
		/** Check to see if convergence has been met.
		 * Arguments:
		 * 		void
		 * Returns:
		 * 		void
		 */

		/** Class Members **/
		mat data;

		int numDat;

		int numCenters;

		int numIter;

		mat centers;

		std::vector<int> newgroup;

		std::vector<int> oldgroup;

		std::vector<int> Rni_init;

};

/**
struct closest{

	closest(int f_indx, float f_dist): indx(f_indx), dist(f_dist){}
	int indx;
	float dist;
};
**/

#endif /* KMEANS_H_ */
