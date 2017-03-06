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

	mean(f_dim[0]*f_dim[1],f_n_models,fill::zeros), noise_var(f_n_models,0), data_dim(3)//, noise_var_model(f_n_models,0)
{
	data = f_data;

	data_dim = f_dim;
}

void PPCA::normalize_data(void){
	std::vector<double> dat_mean(n_var,0);
	std::vector<double> dat_std(n_var, 0);

	for (int i=0; i<n_var; i++){
		for (int j=0; j<n_obs; j++){
			dat_mean[i] += data(i,j);
		}
		dat_mean[i] /= (double)n_obs;
	}


	for (int i=0; i<n_var; i++){
		for (int j=0; j<n_obs; j++){
			dat_std[i] += std::pow(data(i,j) - dat_mean[i], 2);
		}

		dat_std[i] = std::sqrt( dat_std[i] / (double)n_obs);
	}

	for(int i=0; i<n_var; i++){
		for (int j=0; j<n_obs; j++){
			data(i,j) = (data(i,j) - dat_mean[i]) / dat_std[i];
		}
	}

}

int PPCA::get_n_obs(void){

	return n_obs;
}

int PPCA::get_n_var(void){

	return n_var;
}

int PPCA::get_n_components(void){

	return n_components;
}

int PPCA::get_n_models(void){

	return n_models;
}

mat PPCA::get_mean(void){

	return mean;
}

std::vector<double> PPCA::get_noise_variance(void){

	return noise_var;
}

std::vector<int> PPCA::get_data_dim(void){

	return data_dim;
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

	// int status													// TODO: keep track of write success.

	fwrite(Rni.memptr(), sizeof(double),n_obs*n_models, bmpOutput);

	for (int i=0; i<n_models; i++)
		fwrite(&mixfrac[i], sizeof(double), 1, bmpOutput);

	fwrite(mean.memptr(), sizeof(double),n_models*n_var, bmpOutput);

	for (int i=0; i<n_models; i++)
		fwrite(W_mat_vector[i].memptr(), sizeof(double),n_var*n_components, bmpOutput);

	for (int i=0; i<n_models; i++)
		fwrite(S_mat_vector[i].memptr(), sizeof(double),n_var*n_var, bmpOutput);

	for (int i=0; i<n_models; i++){

		fwrite(Minv_mat_vector[i].memptr(), sizeof(double),n_components*n_components, bmpOutput);
	}

	fclose(bmpOutput);

	return 1;
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

			fread(Rni.memptr(), sizeof(double), n_obs*n_models, bmpInput);

			for (int i=0; i<n_models; i++)
				fread(&mixfrac[i], sizeof(double), 1, bmpInput);

			fread(mean.memptr(), sizeof(double),n_models*n_var, bmpInput);

			for (int i=0; i<n_models; i++)
				fread(W_mat_vector[i].memptr(), sizeof(double),n_var*n_components, bmpInput);

			for (int i=0; i<n_models; i++)
				fread(S_mat_vector[i].memptr(), sizeof(double),n_var*n_var, bmpInput);

			for (int i=0; i<n_models; i++)
				fread(Minv_mat_vector[i].memptr(), sizeof(double),n_components*n_components, bmpInput);

			fclose(bmpInput);}

int PPCA_Mixture_EM::write_to_file_mean(char* filename){ //TODO: we need to get the dimensions of the orignial data

	try{

		 mrcParser writeParser = mrcParser(filename,0);

		 float * f_mean_temp = new float[n_var*n_models];

		 for (int i=0;i<mean.n_elem;i++)			// copy from mat to array.
			 f_mean_temp[i] = mean(i);				// armadillo mat is column major so this should work.

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

	for (int n=0; n<n_obs; n++){							// Assign posterior responsibility of mixture i

		for (int i=0; i<n_models; i++)						// for generating data point tn uniformly.
			Rni(n,i) = 1/(double)n_models;
	}
	//normalize_data();
	initialize_helper();									// Initialize other parameters.
}

void PPCA_Mixture_EM::initialize_random(void){

	std::srand(std::time(NULL));							// generate seed.

	for (int n=0; n<n_obs; n++){

		int * temp = new int[n_obs];

		int sum(0);

		for (int i=0; i<n_models; i++){

			temp[i] = std::rand() % 100;					// generate number between 1 and 10
			sum += temp[i];
		}

		for (int i=0; i<n_models; i++)
			Rni(n,i) = temp[i]/(double)sum;

		delete [] temp; temp = NULL;

	}

	//normalize_data();
	initialize_helper();									// Initialize other parameters.
}

void PPCA_Mixture_EM::initialize_kmeans(void){

	kMeans classifier(n_models, 20, data);

	classifier.initCentersRandom();

	std::vector<int> solution;

	solution = classifier.findSolution();

	for (int n=0; n<n_obs; n++)
		Rni(n,solution[n]) = 1;						// Assign total probability to kMeans cluster.
													// Maybe spread probability by gaussian decay around assigned cluster?
	//normalize_data();
	initialize_helper();							// Initialize other parameters.
}

void PPCA_Mixture_EM::initialize_helper(void){

	// The initialize functions above only assign values to Rni. From this, the rest of the parameters are updated.

	for (int i=0; i<n_models; i++){

		mixfrac[i] = 0;															// Initialize model prior probability.

		for (int n=0; n<n_obs; n++)
			mixfrac[i] += Rni(n,i);

		mixfrac[i] /= n_obs;

		colvec temp_mean(n_var, fill::zeros);									// Initialize mean parameter.

		double temp_den = 0;

		for (int n=0; n<n_obs; n++){

			temp_mean += Rni(n,i) * data.col(n);
			temp_den += Rni(n,i);
		}

		mean.col(i) = temp_mean/temp_den;

		mat temp_Si(n_var, n_var, fill::zeros);									// Initialize Si matrices.

		for (int n=0; n<n_obs; n++){

			colvec tempvec = data.col(n) - mean.col(i);
			temp_Si += tempvec * Rni(n,i) * tempvec.t();
		}

		S_mat_vector[i] = temp_Si/((double)mixfrac[i]*n_obs);

		// Update noise variance and W matrices by eigen decomposition of corresponding Si matrix. Iterations after initialization will update by EM rules.
		mat Si_eigvec; 															// Should we designate size?
		vec Si_eigval;

		arma::eig_sym(Si_eigval,Si_eigvec,S_mat_vector[i]);						// Get eigen decomposition.

		double temp_noise_var(0);

		for (int n=0; n< n_var - n_components; n++)								// Initialize noise_var.
			temp_noise_var += Si_eigval(n);										// eig_sym output eigenvalues in ascenting order

		noise_var[i] = temp_noise_var / ((double)(n_var - n_components));		// divide temp_noise_var by d-q

		for (int n=0; n<n_components; n++)										// Initialize W matrices.
			W_mat_vector[i].col(n) = Si_eigvec.col(n_var - n - 1) * std::sqrt((Si_eigval(n_var - n - 1) - noise_var[i]));  // U*sqrt(Aq-sig^2*I)*R

		mat tempM = W_mat_vector[i].t() * W_mat_vector[i];						// Initialize Minv matrices.

		for (int n=0; n<n_components; n++)
			tempM(n,n) += noise_var[i];

		Minv_mat_vector[i] = inv(tempM);
	}

	// All parameters initialized, ready to start iterating.
	has_init = 1;

}

void PPCA_Mixture_EM::update_Rni_all( void){

	double * total = new double[n_obs]();  //doubletotal(0);

	mat log_Ptn_i(n_obs,n_models,fill::zeros);

	for (int i=0; i<n_models; i++){

		mat Cinv = (eye<mat>(n_var,n_var) - W_mat_vector[i] * Minv_mat_vector[i] * W_mat_vector[i].t()) / noise_var[i];

		mat C = eye<mat>(n_var,n_var) * noise_var[i] +  W_mat_vector[i] * W_mat_vector[i].t();

		mat Q_C; mat R_C; mat sign;

		qr(Q_C, R_C, C);											// QR decomposition of C.

		double log_det_C(0);
		for (int k=0; k<n_var; k++)
			log_det_C += std::log(R_C(k,k));

		for (int n=0; n<n_obs; n++){
			Rni(n,i) = calc_log_Ptn_i(n,i,Cinv, log_det_C);
			if (std::isinf(log_Ptn_i(n,i)))
				std:cout << "INF detected log_Prn_i\n";
		}

		double * max = new double[n_obs]();							// get largest log_Ptn_i for each observation
		for (int n=0; n<n_obs; n++)
			max[n] = log_Ptn_i.row(n).max();

		for (int n=0; n<n_obs; n++){								// subtract largest log_Ptn_i from each model for each observation
			for (int i=0; i<n_models; i++){							// This is all for avoiding large numbers. Better numerical stability.
				log_Ptn_i(n,i) -= max[n];
			}
		}
	}

	for (int n=0; n<n_obs; n++){
		for(int i=0; i<n_models; i++){
			Rni(n,i) = std::exp(log_Ptn_i(n,i));
			total[n] += Rni(n,i);
		}
	}

	for (int n=0; n<n_obs; n++){									// Normalize.
		for(int i=0; i<n_models; i++){
			Rni(n,i) *= (mixfrac[n] / total[n]);
		}
	}


	delete [] total; total = NULL;
}

double PPCA_Mixture_EM::calc_Ptn_i(int f_n, int f_i, mat f_Cinv, double f_det_C){			//TODO this function is missing a constant that drops out in the Rni equation. Fix?

	colvec temp = data.col(f_n) - mean.col(f_i);

	double arg = as_scalar(temp.t()* f_Cinv * temp) / 2.0;

	return std::sqrt(1.0/f_det_C) * std::exp(-arg); 	//TODO this is numerically unstable

}

double PPCA_Mixture_EM::calc_log_Ptn_i(int f_n, int f_i, mat f_Cinv, double f_log_det_C){
	colvec temp = data.col(f_n) - mean.col(f_i);

	double arg = as_scalar(temp.t()* f_Cinv * temp) / 2.0;

	double log_Ptn_i = -0.5 * f_log_det_C - arg;

	return log_Ptn_i;
}

void PPCA_Mixture_EM::optimize(int f_max_iter){

	if (!has_init)									// Throw exception if calling this function before initializing parameters.
		throw no_init_exception();

	n_iter = 0;								// reset

	while (n_iter < f_max_iter){

		update_Rni_all();

		for (int i=0; i<n_models; i++){

			mixfrac[i] = 0.0;								// Update model prior probability.

			for (int n=0; n<n_obs; n++)
				mixfrac[i] += Rni(n,i);

			mixfrac[i] /= (double)n_obs;
		}

		for (int i=0; i<n_models; i++){

			colvec temp_mean(n_var,fill::zeros);			// Update mean parameter.
			double temp_den = 0;

			for (int n=0; n<n_obs; n++){

				temp_mean += Rni(n,i) * data.col(n);
				temp_den += Rni(n,i);
			}

			mean.col(i) = temp_mean/temp_den;
		}

		for (int i=0; i<n_models; i++){						// Update W mat, noise_var, S mat, and Minv mat last

			mat tempSW_mat(n_var, n_components, fill::zeros);			//Calculate S*W.
			/**
			 * The below lines are meant to calculate S*W efficiently, but we need to divide by
			 * mixfrac[i]*n_obs. Unfortunately, mixfrac is updated first. We need to save a copy
			 * of the old mixfrac value before update.
			for (int n=0; n<n_obs; n++){

				colvec tempvec = data.col(n) - mean(i);
				tempSW_mat += Rni(n,i) * tempvec * (tempvec.t() * W_mat_vector[i]);
			}
			**/
			tempSW_mat =  S_mat_vector[i] *W_mat_vector[i];


			mat tempInv_mat = inv(noise_var[i] * eye(n_components, n_components) + Minv_mat_vector[i] * W_mat_vector[i].t() * tempSW_mat);

			W_mat_vector[i] = tempSW_mat * tempInv_mat;					// Update W mat.

			noise_var[i] = 1.0/(double)n_var * trace(S_mat_vector[i] - tempSW_mat * Minv_mat_vector[i] * W_mat_vector[i].t());

			mat temp_Si(n_var,n_var, fill::zeros);									// Update S matrices.

			for (int n=0; n<n_obs; n++){

				colvec tempvec = data.col(n) - mean.col(i);
				temp_Si += Rni(n,i)*tempvec*tempvec.t();
			}

			S_mat_vector[i] = temp_Si/(mixfrac[i] * n_obs);

			mat tempM = eye<mat>(n_components, n_components);

			Minv_mat_vector[i] = inv(noise_var[i] * eye<mat>(n_components, n_components) + W_mat_vector[i].t() * W_mat_vector[i]); // Update Minv matrices.

		}

		std::cout << "Iteration " << n_iter + 1 << " completed.\n";
		this->n_iter++;

		if (!n_iter % 10 ){
			std::stringstream ss;
			ss << "params_iter_" << n_iter << ".bmp";
			std::string outfile_s = ss.str();
			char * outFile_c = new char[outfile_s.length() + 1];	// TODO: need to get rid of using *char as input.
			strcpy(outFile_c, outfile_s.c_str());
			write_to_file_params(outFile_c);
			delete [] outFile_c; outFile_c = NULL;
		}

		std::stringstream ss;
		ss << "mean_iter_" << n_iter << ".mrc";
		std::string outfile_s = ss.str();
		char * outFile_c = new char[outfile_s.length() + 1];
		strcpy(outFile_c, outfile_s.c_str());
		write_to_file_mean(outFile_c);
		delete [] outFile_c; outFile_c = NULL;

	}

}



