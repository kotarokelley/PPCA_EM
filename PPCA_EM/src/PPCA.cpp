//============================================================================
// Name        : PPCA.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include "PPCA.h"

using namespace arma;

PPCA_Mixture_EM::PPCA_Mixture_EM(fmat f_data, int* f_dim):

	n_components(f_dim[1]),		// Set to number of data points so keep all components.

	components_(1,fmat(f_dim[0],f_dim[1])), explained_variance_ratio(1,frowvec(f_dim[1])), mean(1,frowvec(f_dim[0])),

	noise_variance(1,0), noise_variance_model(1,0), n_models(1), W_mat_vector(1,0), Si_mat_vector(1,0)

{
	mixfrac = new float[1];
	Rni = new float[1];
	Rn = new float[1];

	data = f_data;

}

PPCA_Mixture_EM::PPCA_Mixture_EM(fmat f_data, int* f_dim, int f_n_components):
	n_components(f_n_components),			// Keep only the largest components.
	components_(1,fmat(f_dim[0],f_n_components)), explained_variance_ratio(1,frowvec(f_n_components)), mean(1,frowvec(f_dim[0])),
	noise_variance(1,0), noise_variance_model(1,0), n_models(1), W_mat_vector(1,0), Si_mat_vector(1,0)
{

	mixfrac = new float[1];
	Rni = new float[1];
	Rn = new float[1];

	data = f_data;

}

PPCA_Mixture_EM::PPCA_Mixture_EM(fmat f_data, int* f_dim, int f_n_components, int f_n_models):
		n_components(f_n_components),		// Keep only the largest components.
		components_(f_n_models,fmat(f_dim[0],f_n_components)), explained_variance_ratio(f_n_models,frowvec(f_n_components)), mean(f_n_models,frowvec(f_dim[0])),
		noise_variance(1,0), noise_variance_model(1,0), n_models(f_n_models), W_mat_vector(1,0), Si_mat_vector(1,0){

	mixfrac = new float[1];
	Rni = new float[1];
	Rn = new float[1];

	data = f_data;

}

void PPCA_Mixture_EM::fit(void){

}

char * PPCA_Mixture_EM::get_params(void){
	return "eoiheoih";
}

int PPCA_Mixture_EM::get_n_components(void){
	return this->n_components;
}

std::vector<fmat> PPCA_Mixture_EM::get_components_(void){
	return this->components_;
}

std::vector<frowvec> PPCA_Mixture_EM::get_explained_variance_ratio(void){
	return this->explained_variance_ratio;
}
std::vector<frowvec> PPCA_Mixture_EM::get_mean(void){
	return this->mean;
}

std::vector<frowvec> PPCA_Mixture_EM::get_noise_variance(void){
	return noise_variance;
}

std::vector<frowvec> PPCA_Mixture_EM::get_noise_variance_moedl(void){
	return noise_variance_model;
}

