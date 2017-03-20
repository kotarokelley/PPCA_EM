/*
 * kMeans.cpp
 *
 *  Created on: Nov 13, 2016
 *      Author: kotarokelley
 */

#include "kMeans.h"

using namespace arma;

kMeans::kMeans(int f_numCenters, int f_numIter, mat f_data):

	data(f_data), numDat(f_data.n_cols), numCenters(f_numCenters), numIter(f_numIter), centers(f_data.n_rows,f_numCenters),

	newgroup(f_data.n_cols,-2), oldgroup(f_data.n_cols,-1), Rni_init(f_numCenters)
{

}

void kMeans::initCentersRandom(void){

	std::srand(std::time(NULL));										// generate seed.
	int * indices = new int[numDat]();
	int counter = 0;

	while (counter<numCenters){
		int indx = std::rand() % numDat;							// randomly choose indices for centers. draw from 0-(data.n_cols-1)
		if (!indices[indx]){
			indices[indx] = 1;
			centers.col(counter) = data.col(indx);
			counter++;
		}
	}

	delete [] indices; indices = NULL;
}

std::vector<int> kMeans::findSolution(void){

	int counter = 0;

	while(counter <numIter){//!hasConverged() && counter < numIter){			// keep looping until return conditions are met.

		for (int i=0; i<numDat; i++)
			oldgroup[i] = newgroup[i];						// copy newgroup to oldgroup before next iteration of center assignment for each data point

		for (int i=0; i<numDat; i++)
			newgroup[i] = findClosest(i);// Find closest center for each data point.

		recalc_Centers();
		counter++;
	}

	return newgroup;												// Converged or maximum iterations met.
}

mat kMeans::getCenters(void){
	return centers;
}

int kMeans::findClosest(int in_indx){

	double dist = std::numeric_limits<double>::max();						//set to max double
	double temp;
	int out_indx;

	for (int i=0; i<numCenters; i++){
		temp = norm(centers.col(i) - data.col(in_indx),2);   // calculate distance to each center

		if (temp < dist){
			dist = temp;
			out_indx = i;
		}
	}

	return out_indx;
}

void kMeans::recalc_Centers(void){

	int * counts = new int[numCenters]();

	mat tempCenters(centers.n_rows,centers.n_cols,fill::zeros);

	for (int i=0; i<numDat; i++){

		tempCenters.col(newgroup[i]) += data.col(i);
		counts[newgroup[i]] ++;
	}

	for (int i=0; i<numCenters; i++){
		if(counts[i]>=1)
			centers.col(i) = tempCenters.col(i)/ (double)counts[i];
		else
			std::cout << i << " ERROR: kmeans counts for center: " << i << " is 0" << "\n"; // TODO: throw exception
	}
	delete [] counts; counts = NULL;
}

bool kMeans::hasConverged(void){
	bool converged = true;
	for (int i=0; i<numDat; i++){
		if (newgroup[i] != oldgroup[i])				// compare indices from the new and old centers for each data point
			converged = false;
	}
	return converged;
}



