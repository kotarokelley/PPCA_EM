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
	int * indices = new int[this->numDat];
	int counter = 0;

	while (counter<this->numCenters){
		int indx = std::rand() % this->numDat;							// randomly choose indices for centers. draw from 0-(data.n_cols-1)
		if (!indices[indx]){
			indices[indx] = 1;
			this->centers.col(counter) = this->data.col(indx);
			counter++;
		}
	}

	delete [] indices; indices = NULL;
}

std::vector<int> kMeans::findSolution(void){
	int counter = 0;
	//float maxDist = std::numeric_limits<float>::max();					// set to max float


	while(!this->hasConverged() && counter < this->numIter){			// keep looping until return conditions are met.
		for (int i=0; i<this->numDat; i++){
			this->oldgroup[i] = this->newgroup[i];						// copy newgroup to oldgroup before next iteration of center assignment for each data point
		}

		for (int i=0; i<this->numDat; i++){
			int closest_indx = this->findClosest(i);					// Find closest center for each data point.
			this->newgroup[i] = closest_indx;
		}

		this->recalc_Centers();
		counter++;
	}

	return this->newgroup;												// Converged or maximum iterations met.

}

int kMeans::findClosest(int in_indx){
	double dist = std::numeric_limits<float>::max();						//set to max float
	double temp;
	int out_indx;
	for (int i=0; i<this->numCenters; i++){
		temp = norm_dot(this->centers.col(i),this->data.col(in_indx));   // calculate distance to each center
		if (temp < dist){
			dist = temp;
			out_indx = i;
		}
	}

	return out_indx;

}

void kMeans::recalc_Centers(void){
	mat temp(this->data);
	for (int i=0; i<this->numCenters; i++){

	}
}

bool kMeans::hasConverged(void){
	bool converged = true;
	for (int i=0; i<this->numDat; i++){
		if (this->newgroup[i] != this->oldgroup[i])				// compare indices from the new and old centers for each data point
			converged = false;
	}
	return converged;
}



