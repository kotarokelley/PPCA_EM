//============================================================================
// Name        : PPCA.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include "PPCA.h"

using namespace arma;

PPCA::PPCA(mat f_data, int* f_dim, int f_n_components, int f_n_models):

	n_obs(f_dim[1]), n_var(f_dim[0]),

	n_components(f_n_components), n_models(f_n_models),	// Keep only the largest components.

	components_(f_n_models,mat(f_dim[0],f_n_components)), explained_variance_ratio(f_n_models,rowvec(f_n_components)),

	mean(f_dim[0],f_n_components,fill::zeros), noise_var(f_n_models,0)//, noise_var_model(f_n_models,0)
{
	data = f_data;
}


PPCA_Mixture_EM::PPCA_Mixture_EM(mat f_data, int* f_dim, int f_n_components, int f_n_models):

		PPCA(f_data, f_dim, f_n_components, f_n_models),

		mixfrac(f_n_models,0), n_iter(0), Rni(f_dim[1],f_n_models,fill::zeros),

		W_mat_vector(f_n_models,mat(f_dim[0],f_n_components,fill::zeros)), Si_mat_vector(f_n_models,mat(f_dim[0],f_dim[0],fill::zeros)),

		Minv_mat_vector(f_n_models,mat(f_n_components,f_n_components,fill::zeros)), has_init(0)
{

}

void PPCA_Mixture_EM::fit(void){

}

std::tuple<int, int, int, int> PPCA_Mixture_EM::get_params(void){

	return std::make_tuple(this->n_var,  this->n_obs, this->n_components, this->n_models);
}

int PPCA_Mixture_EM::get_n_components(void){
	return this->n_components;
}

std::vector<mat> PPCA_Mixture_EM::get_components_(void){
	return this->components_;
}

mat PPCA_Mixture_EM::get_mean(void){
	return this->mean;
}

std::vector<double> PPCA_Mixture_EM::get_noise_variance(void){
	return this->noise_var;
}


void PPCA_Mixture_EM::initialize_uniform(void){
	for (int n=0; n<this->n_obs; n++){					// Assign posterior responsibility of mixture i
		for (int i=0; i<this->n_models; i++){					// for generating data point tn uniformly.
			Rni(n,i) = 1/(double)this->n_models;
		}
	}

	this->initialize_helper();								// Initialize other parameters.

}

void PPCA_Mixture_EM::initialize_random(void){
	std::srand(std::time(NULL));							// generate seed.


	for (int n=0; n<this->n_obs; n++){

		int * temp = new int[this->n_obs];
		int sum(0);

		for (int i=0; i<this->n_models; i++){
			temp[i] = std::rand() % 100;					// generate number between 1 and 10
			sum += temp[i];
		}

		for (int i=0; i<this->n_models; i++)
			Rni(n,i) = temp[i]/(double)sum;

		delete [] temp; temp = NULL;

	}

	this->initialize_helper();								// Initialize other parameters.
}

void PPCA_Mixture_EM::initialize_kmeans(void){

	kMeans classifier(this->n_models,20,this->data);
	classifier.initCentersRandom();
	std::vector<int> solution;
	solution = classifier.findSolution();

	for (int n=0; n<this->n_obs; n++){
		this->Rni(n,solution[n]) = 1;				// Assign total probability to kMeans cluster.
													// Maybe spread probability by gaussian decay around assigned cluster?
	}

	this->initialize_helper();						// Initialize other parameters.

}

void PPCA_Mixture_EM::initialize_helper(void){
	// The initialize functions above only assign values to Rni. From this, the rest of the parameters are updated.
	for (int i=0; i<this->n_models; i++){			// TODO: Change the limit to this->n_models

		this->mixfrac[i] = 0;								// Initialize model prior probability.
		for (int n=0; n<this->n_obs; n++){
			this->mixfrac[i] += this->Rni(n,i);
		}
		this->mixfrac[i] /= this->n_obs;

		colvec temp_mean(this->n_var,fill::zeros);			// Initialize mean parameter.
		double temp_den = 0;
		for (int n=0; n<this->n_obs; n++){
			temp_mean += this->Rni(n,i)*this->data.col(n);
			temp_den += this->Rni(n,i);
		}
		this->mean.col(i) = temp_mean/temp_den;

		mat temp_Si(this->n_var,this->n_obs,fill::zeros);									// Initialize Si matrices.
		for (int n=0; n<this->n_obs; n++){
			colvec tempvec = this->data.col(n) - this->mean.col(i);
			temp_Si = this->Rni(n,i)*tempvec*tempvec.t();
			//temp_Si += this->Rni(n,i)*(this->data.col(n) - this->mean.col(i))*(this->data.col(n) - this->mean.col(i)).t();
		}
		this->Si_mat_vector[i] = temp_Si/(this->mixfrac[i]*this->n_obs);

		// Update noise variance and W matrices by eigen decomposition of corresponding Si matrix. Iterations after initialization will update by EM rules.
		mat Si_eigvec; 				// Function eig_gen takes in a col<cx_double>, mat<cx_double>, mat<double>
		vec Si_eigval;
		arma::eig_sym(Si_eigval,Si_eigvec,this->Si_mat_vector[i]);		// Get eigen decomposition.
		double temp_noise_var(0);
		for (int n=this->n_components; n<this->n_var; n++){									// Initialize noise_var.
			temp_noise_var += Si_eigval(n);
		}
		this->noise_var[i] = temp_noise_var / (this->n_var - this->n_components);			// divide temp_noise_var by d-q

		mat temp_diag_mat(this->n_components,this->n_components,fill::zeros);				// Initialize W matrices.
		for (int n=0; n<this->n_components;n++){
			temp_diag_mat(n,n) = std::sqrt((Si_eigval(n) - this->noise_var[i]));			//sqrt(Aq-sig^2*I)
		}
		this->W_mat_vector[i] = (Si_eigvec.cols(0,this->n_components-1))*temp_diag_mat;		// U*sqrt(Aq-sig^2*I)*R

		mat tempM = eye<mat>(this->n_components,this->n_components);
		this->Minv_mat_vector[i] = this->noise_var[i]*(inv(this->noise_var[i]*eye<mat>(this->n_components, 		// Initialize Minv matrices.
				this->n_components) + this->W_mat_vector[i].t()*this->W_mat_vector[i]));

		// All parameters initialized, ready to start iterating.
	}
	//for (int i=0; i<this->n_models; i++){
		//std::cout<<this->noise_var[i]<<"\n";
	//}
	this->has_init = 1;

}

void PPCA_Mixture_EM::update_Rni(int n){
	double total(0);
	for (int i=0; i<this->n_models; i++){
		double temp_Ptn_i = this->calc_Ptn_i(n,i);
		total += this->Rni(n,i);
		this->Rni(n,i) = temp_Ptn_i * this->mixfrac[i];
		}
	for (int i=0; i<this->n_models; i++){
		Rni(n,i) /= total;
	}
}

double PPCA_Mixture_EM::calc_Ptn_i(int n, int i){			//TODO this function is missing a constant that drops out in the Rni equation. Fix?
	mat Cinv = eye<mat>(this->n_var,this->n_var) * std::pow(this->noise_var[i],-2.0) -
			this->W_mat_vector[i]*this->Minv_mat_vector[i]*this->W_mat_vector[i].t() * std::pow(this->noise_var[i],-4.0);

	colvec temp = this->data.col(n) - this->mean.col(i);
	double arg = as_scalar(temp.t()* Cinv* temp);
	double deter = det(Cinv);

	return std::pow(deter,-0.5) * std::exp(arg);

}



void PPCA_Mixture_EM::optimize(int f_max_iter){

	if (!this->has_init)									// Throw exception if calling this function before initializing parameters.
		throw no_init_exception();

	while (!(this->n_iter <= f_max_iter)){

		for (int n=0; n<this->n_obs; n++){						// Update Rni first.
			this->update_Rni(n);
		}

		for (int i=0; i<this->n_models; i++){
			this->mixfrac[i] = 0.0;								// Update model prior probability.
			for (int n=0; n<this->n_obs; n++){
				this->mixfrac[i] += this->Rni(n,i);
			}
			this->mixfrac[i] /= this->n_obs;
		}

		for (int i=0; i<this->n_models; i++){
			colvec temp_mean(this->n_var,fill::zeros);			// Update mean parameter.
			double temp_den = 0;
			for (int n=0; n<this->n_obs; n++){
				temp_mean += this->Rni(n,i)*this->data.col(n);
				temp_den += this->Rni(n,i);
			}
			this->mean.col(i) = temp_mean/temp_den;
		}


		for (int i=0; i<this->n_models; i++){					// Update W mat, noise_var, S mat, and Minv mat last

			mat tempSW_mat(this->n_components,this->n_var,fill::zeros);			//Calculate S*W.
			for (int n=0; n<this->n_obs; n++){
				colvec tempvec = this->data.col(n) - this->mean(i);
				tempSW_mat += this->Rni(n,i) * tempvec * (tempvec.t() * this->W_mat_vector[i]);
			}
			tempSW_mat /= this->mixfrac[i] * this->n_obs;

			mat tempInv_mat = inv( this->noise_var[i] * eye(this->n_components,this->n_components) +
					this->Minv_mat_vector[i] * this->W_mat_vector[i].t() * tempSW_mat);

			this->W_mat_vector[i] = tempSW_mat * tempInv_mat;					// Update W mat.

			this->noise_var[i] = 1.0/this->n_var * trace(this->Si_mat_vector[i] // Update noise_var.
									- tempSW_mat * this->Minv_mat_vector[i] * this->W_mat_vector[i].t());


			mat temp_Si(this->n_var,this->n_obs,fill::zeros);									// Update S matrices.
			for (int n=0; n<this->n_obs; n++){
				colvec tempvec = this->data.col(n) - this->mean.col(i);
				temp_Si = this->Rni(n,i)*tempvec*tempvec.t();
			}
			this->Si_mat_vector[i] = temp_Si/(this->mixfrac[i]*this->n_obs);

			mat tempM = eye<mat>(this->n_components,this->n_components);
			this->Minv_mat_vector[i] = this->noise_var[i]*(inv(this->noise_var[i]*eye<mat>(this->n_components, 		// Initialize Minv matrices.
					this->n_components) + this->W_mat_vector[i].t()*this->W_mat_vector[i]));


		}



		this->n_iter++;
	}

}



