//============================================================================
// Name        : PPCA.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include "PPCA.h"

using namespace arma;

PPCA::PPCA(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models):

	n_obs(f_dim[2]), n_var(f_dim[0]*f_dim[1]),

	n_components(f_n_components), n_models(f_n_models),	// Keep only the largest components.

	//components_(f_n_models,mat(f_dim[0],f_n_components)), explained_variance_ratio(f_n_models,rowvec(f_n_components)),

	mean(f_dim[0]*f_dim[1],f_n_components,fill::zeros), noise_var(f_n_models,0), data_dim(3)//, noise_var_model(f_n_models,0)
{
	data = f_data;

	data_dim = f_dim;

	//data_dim[0] = f_dim[0];			// This is kind of clunky. Should we do data_dim.push_back()?

	//data_dim[1] = f_dim[1];			// can we just say data_dim = f_dim?

	//data_dim[2] = f_dim[2];

}


PPCA_Mixture_EM::PPCA_Mixture_EM(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models):

		PPCA(f_data, f_dim, f_n_components, f_n_models),

		mixfrac(f_n_models,0), n_iter(0), Rni(f_dim[2],f_n_models,fill::zeros),

		W_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros)), S_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_dim[0]*f_dim[1],fill::zeros)),

		Minv_mat_vector(f_n_models,mat(f_n_components,f_n_components,fill::zeros)), has_init(0)
{

}

PPCA_Mixture_EM::PPCA_Mixture_EM(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, char* prev_run_filename):

		PPCA(f_data, f_dim, f_n_components, f_n_models),

		mixfrac(f_n_models,0), n_iter(0), Rni(f_dim[2],f_n_models,fill::zeros),

		W_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros)), S_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_dim[0]*f_dim[1],fill::zeros)),

		Minv_mat_vector(f_n_models,mat(f_n_components,f_n_components,fill::zeros)), has_init(0)
{

	parse_prev_params(prev_run_filename);			//TODO this is doing redundant work to the initalization list.

	has_init = 1;
}

int PPCA_Mixture_EM::get_n_components(void){

	return n_components;
}


mat PPCA_Mixture_EM::get_mean(void){

	return mean;
}

std::vector<double> PPCA_Mixture_EM::get_noise_variance(void){

	return noise_var;
}

std::vector<int> PPCA_Mixture_EM::get_data_dim(void){

	return data_dim;
}

int PPCA_Mixture_EM::write_to_file_params(char* filename){

	FILE * bmpOutput = fopen(filename, "wb");		// OPEN FILE

	if(bmpOutput == NULL){

		printf("Could not open file: %s", filename);
		exit (1);
	}

	fseek(bmpOutput, 0, SEEK_SET);		    // SET POINTER TO BEGINNING OF FILE


	fwrite(&data_dim[0], sizeof(int),1,bmpOutput);		// xpix

	fwrite(&data_dim[1], sizeof(int),1,bmpOutput);		// ypix

	fwrite(&data_dim[2], sizeof(int),1,bmpOutput);			// n particles

	fwrite(&n_models, sizeof(int),1,bmpOutput);			// n models

	fwrite(&n_components, sizeof(int),1,bmpOutput);		// n components

	double * Rni_temp = new double[n_obs*n_models];		// Copy Rni. Not sure if this can be done another way.

	for (int i=0; i<n_obs*n_models; i++){

		Rni_temp[i] = Rni(i);							// armadillo is column major.
	}

	// int status													// TODO: keep track of write success.

	fwrite(Rni_temp, sizeof(double),n_obs*n_models, bmpOutput);

	delete [] Rni_temp; Rni_temp = NULL;

	double * mixfrac_temp = new double[n_models];

	for (int i=0; i<n_models; i++){

		mixfrac_temp[i] = mixfrac[i];
	}

	fwrite(mixfrac_temp, sizeof(double), n_models, bmpOutput);

	delete [] mixfrac_temp; mixfrac_temp = NULL;
	/**
	double * W_mat_vector_temp = new double[n_models*n_components*n_var];

	for (int i=0; i<n_models; i++){
		for (int j=0; j<n_components*n_var; j++){
			W_mat_vector_temp[i] = W_mat_vector[i](j);
		}
	}

	fwrite(W_mat_vector_temp, sizeof(double),n_models*n_components*n_var, bmpOutput);
	delete [] W_mat_vector_temp; W_mat_vector_temp = NULL;

	double * S_mat_vector_temp = new double[n_var*n_var*n_models];

	for (int i=0; i<n_models; i++){
		for (int j=0; j<n_var*n_var; j++){
			S_mat_vector_temp[i] = S_mat_vector[i](j);
		}
	}

	fwrite(S_mat_vector_temp, sizeof(double),n_models*n_var*n_var, bmpOutput);
	delete [] S_mat_vector_temp; S_mat_vector_temp = NULL;


	double * Minv_mat_vector_temp = new double[n_components*n_components];
	**/

	fwrite(mean.memptr(), sizeof(double),n_models*n_var, bmpOutput);

	for (int i=0; i<n_models; i++){

		fwrite(W_mat_vector[i].memptr(), sizeof(double),n_var*n_components, bmpOutput);
	}

	for (int i=0; i<n_models; i++){

		fwrite(S_mat_vector[i].memptr(), sizeof(double),n_var*n_var, bmpOutput);
	}

	for (int i=0; i<n_models; i++){

		fwrite(Minv_mat_vector[i].memptr(), sizeof(double),n_components*n_components, bmpOutput);
	}

	fclose(bmpOutput);

}

void PPCA_Mixture_EM::parse_prev_params(char* prev_run_filename){

			FILE * bmpInput = fopen(prev_run_filename, "rb");		// OPEN FILE
			if(bmpInput == NULL){
				printf("Could not open file: %s", prev_run_filename);
				exit (1);
			}

			fseek(bmpInput, 0, SEEK_END);			// SET POINTER TO END OF FILE

			unsigned long fileSize = ftell (bmpInput);			// GET TOTAL FILE SIZE IN BYTES

			fseek(bmpInput, 0, SEEK_SET);		    // SET POINTER TO BEGINNING OF FILE

			fread(&data_dim[0], sizeof(int),1,bmpInput);		// xpix

			fread(&data_dim[1], sizeof(int),1,bmpInput);		// ypix

			fread(&data_dim[2], sizeof(int),1,bmpInput);		// ypix

			n_obs = data_dim[2];

			n_var = data_dim[0]*data_dim[1];

			fread(&n_models, sizeof(int),1,bmpInput);			// n models

			fread(&n_components, sizeof(int),1,bmpInput);		// n components

			fread(mean.memptr(), sizeof(double),n_models*n_var, bmpInput);

			for (int i=0; i<n_models; i++){

				fread(W_mat_vector[i].memptr(), sizeof(double),n_var*n_components, bmpInput);
			}

			for (int i=0; i<n_models; i++){

				fread(S_mat_vector[i].memptr(), sizeof(double),n_var*n_var, bmpInput);
			}

			for (int i=0; i<n_models; i++){

				fread(Minv_mat_vector[i].memptr(), sizeof(double),n_components*n_components, bmpInput);
			}

			fclose(bmpInput);}

int PPCA_Mixture_EM::write_to_file_mean(char* filename){ //TODO: we need to get the dimensions of the orignial data

	try{

		 mrcParser writeParser = mrcParser(filename);

		 float * f_mean_temp = new float[n_var*n_models];

		 for (int i=0;i<mean.n_elem;i++){			// copy from mat to array.

			 f_mean_temp[i] = mean(i);				// armadillo mat is column major so this should work.
		 }

		 std::vector<int> out_dim(3);

		 out_dim[0] = data_dim[0]; out_dim[1] = data_dim[1]; out_dim[2] = n_models;

		 writeParser.writeData(f_mean_temp, out_dim);

		 delete [] f_mean_temp, f_mean_temp = NULL;

		 return 1;
	}
	catch (...){				// TODO: make more specific.

		return 0;
	}

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

	initialize_helper();								// Initialize other parameters.
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

	initialize_helper();						// Initialize other parameters.

}

void PPCA_Mixture_EM::initialize_helper(void){
	// The initialize functions above only assign values to Rni. From this, the rest of the parameters are updated.
	for (int i=0; i<n_models; i++){

		mixfrac[i] = 0;								// Initialize model prior probability.

		for (int n=0; n<n_obs; n++){

			mixfrac[i] += Rni(n,i);
		}

		mixfrac[i] /= n_obs;

		colvec temp_mean(n_var, fill::zeros);			// Initialize mean parameter.

		double temp_den = 0;

		for (int n=0; n<n_obs; n++){

			temp_mean += Rni(n,i) * data.col(n);

			temp_den += Rni(n,i);
		}

		this->mean.col(i) = temp_mean/temp_den;

		mat temp_Si(n_var, n_var, fill::zeros);									// Initialize Si matrices.

		for (int n=0; n<n_obs; n++){

			colvec tempvec = data.col(n) - mean.col(i);

			temp_Si += Rni(n,i)*tempvec*tempvec.t();
		}

		S_mat_vector[i] = temp_Si/((double)mixfrac[i]*n_obs);

		// Update noise variance and W matrices by eigen decomposition of corresponding Si matrix. Iterations after initialization will update by EM rules.
		mat Si_eigvec; 				// Function eig_gen takes in a col<cx_double>, mat<cx_double>, mat<double>

		vec Si_eigval;

		arma::eig_sym(Si_eigval,Si_eigvec,S_mat_vector[i]);		// Get eigen decomposition.

		double temp_noise_var(0);

		for (int n=n_components; n<n_var; n++){									// Initialize noise_var.

			temp_noise_var += Si_eigval(n);
		}

		noise_var[i] = temp_noise_var / (n_var - n_components);			// divide temp_noise_var by d-q

		mat temp_diag_mat(n_components, n_components, fill::zeros);				// Initialize W matrices.

		for (int n=0; n<n_components; n++){

			temp_diag_mat(n,n) = std::sqrt((Si_eigval(n) - noise_var[i]));			//sqrt(Aq-sig^2*I)
		}

		this->W_mat_vector[i] = (Si_eigvec.cols(0, n_components-1))*temp_diag_mat;		// U*sqrt(Aq-sig^2*I)*R

		mat tempM = eye<mat>(n_components, n_components);

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
															//TODO: This function calculates Cinv for each particle for each model. Very inefficient.
double PPCA_Mixture_EM::calc_Ptn_i(int n, int i){			//TODO this function is missing a constant that drops out in the Rni equation. Fix?
	mat Cinv = eye<mat>(this->n_var,this->n_var) * std::pow(this->noise_var[i],-2.0) -
			this->W_mat_vector[i]*this->Minv_mat_vector[i]*this->W_mat_vector[i].t() * std::pow(this->noise_var[i],-4.0);

	colvec temp = this->data.col(n) - this->mean.col(i);
	double arg = as_scalar(temp.t()* Cinv* temp);
	double deter = det(Cinv);

	return std::pow(deter,-0.5) * std::exp(arg);

}

void PPCA_Mixture_EM::optimize(int f_max_iter){

	if (!has_init)									// Throw exception if calling this function before initializing parameters.
		throw no_init_exception();

	n_iter = 0;								// reset

	while (n_iter <= f_max_iter){

		for (int n=0; n<n_obs; n++){						// Update Rni first.

			update_Rni(n);
		}

		for (int i=0; i<n_models; i++){

			mixfrac[i] = 0.0;								// Update model prior probability.

			for (int n=0; n<n_obs; n++){

				this->mixfrac[i] += Rni(n,i);
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

			mat tempSW_mat(this->n_var,this->n_components,fill::zeros);			//Calculate S*W.
			for (int n=0; n<this->n_obs; n++){
				colvec tempvec = this->data.col(n) - this->mean(i);
				tempSW_mat += this->Rni(n,i) * tempvec * (tempvec.t() * this->W_mat_vector[i]);
			}
			tempSW_mat /= this->mixfrac[i] * this->n_obs;

			mat tempInv_mat = inv( this->noise_var[i] * eye(this->n_components,this->n_components) +
					this->Minv_mat_vector[i] * this->W_mat_vector[i].t() * tempSW_mat);

			this->W_mat_vector[i] = tempSW_mat * tempInv_mat;					// Update W mat.

			this->noise_var[i] = 1.0/this->n_var * trace(this->S_mat_vector[i] // Update noise_var.
									- tempSW_mat * this->Minv_mat_vector[i] * this->W_mat_vector[i].t());


			mat temp_Si(this->n_var,this->n_obs,fill::zeros);									// Update S matrices.
			for (int n=0; n<this->n_obs; n++){
				colvec tempvec = this->data.col(n) - this->mean.col(i);
				temp_Si = this->Rni(n,i)*tempvec*tempvec.t();
			}
			this->S_mat_vector[i] = temp_Si/(this->mixfrac[i]*this->n_obs);

			mat tempM = eye<mat>(this->n_components,this->n_components);
			this->Minv_mat_vector[i] = this->noise_var[i]*(inv(this->noise_var[i]*eye<mat>(this->n_components, 		// Update Minv matrices.
					this->n_components) + this->W_mat_vector[i].t()*this->W_mat_vector[i]));


		}



		this->n_iter++;
	}

}



