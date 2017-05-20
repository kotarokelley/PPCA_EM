//============================================================================
// Name        : PPCA.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include "PPCA.h"

#define DEBUG

using namespace arma;

PPCA::PPCA(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models):

	n_obs(f_dim[2]), n_var(f_dim[0]*f_dim[1]),

	n_components(f_n_components), n_models(f_n_models),	// Keep only the largest components.

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

		Minv_mat_vector(f_n_models,mat(f_n_components,f_n_components,fill::zeros)), log_likelihood(LARGE_NUMBER), has_init(0)
{

}

PPCA_Mixture_EM::PPCA_Mixture_EM(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, char* prev_run_filename):

		PPCA(f_data, f_dim, f_n_components, f_n_models),

		mixfrac(f_n_models,0), n_iter(0), Rni(f_dim[2],f_n_models,fill::zeros),

		W_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros)), S_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_dim[0]*f_dim[1],fill::zeros)),

		Minv_mat_vector(f_n_models,mat(f_n_components,f_n_components,fill::zeros)), log_likelihood(LARGE_NUMBER), has_init(0)
{

	parse_prev_params(prev_run_filename);			//TODO this is doing redundant work to the initalization list.

	has_init = 1;
}

int PPCA_Mixture_EM::write_to_file_params(int n_iter){

	std::stringstream ss;

	ss << "params_iter_" << n_iter << ".bmp";

	std::string outfile_s = ss.str();

	char * outFile_c = new char[outfile_s.length() + 1];

	strcpy(outFile_c, outfile_s.c_str());

	FILE * bmpOutput = fopen(outFile_c, "wb");		// OPEN FILE

	if(bmpOutput == NULL){
		std::cout << "Error function write_to_file_params: Could not open file: " << outFile_c << "\n";
		std::cout << "Exiting program.\n";
		exit (1);
	}

	fseek(bmpOutput, 0, SEEK_SET);		    // SET POINTER TO BEGINNING OF FILE

	fwrite(&data_dim[0], sizeof(int),1,bmpOutput);		// xpix

	fwrite(&data_dim[1], sizeof(int),1,bmpOutput);		// ypix

	fwrite(&data_dim[2], sizeof(int),1,bmpOutput);			// n particles

	fwrite(&n_models, sizeof(int),1,bmpOutput);			// n models

	fwrite(&n_components, sizeof(int),1,bmpOutput);		// n components

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
		std::cout << "Error function parse_prev_params: Could not open file: " << prev_run_filename << "\n";
		std::cout << "Exiting program.\n";
		exit (1);
	}

	fseek(bmpInput, 0, SEEK_END);			// SET POINTER TO END OF FILE

	//unsigned long fileSize = ftell (bmpInput);			// GET TOTAL FILE SIZE IN BYTES

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

	fclose(bmpInput);
}

int PPCA_Mixture_EM::write_to_file_mean(int n_iter){

	std::stringstream ss;

	ss << "mean_iter_" << n_iter << ".mrc";

	std::string outfile_s = ss.str();

	char * outFile_c = new char[outfile_s.length() + 1];

	strcpy(outFile_c, outfile_s.c_str());

	try{
		 mrcParser writeParser = mrcParser(outFile_c,0);

		 float * f_mean_temp = new float[n_var*n_models];

		 for (int i=0;i<mean.n_elem;i++)			// copy from mat to array.
			 f_mean_temp[i] = mean(i);				// armadillo mat is column major so this should work.

		 std::vector<int> out_dim(3);

		 out_dim[0] = data_dim[0]; out_dim[1] = data_dim[1]; out_dim[2] = n_models;

		 writeParser.writeData(f_mean_temp, out_dim);

		 delete [] f_mean_temp, f_mean_temp = NULL;
		 delete [] outFile_c; outFile_c = NULL;

		 return 1;
	}
	catch (...){				// TODO: make more specific.

		delete [] outFile_c; outFile_c = NULL;
		return 0;
	}

}


void PPCA_Mixture_EM::initialize_random(void){

	std::srand(std::time(NULL));							// generate seed.

	for (int n=0; n<n_obs; n++){

		double * temp = new double[n_models];

		double sum(0);

		for (int i=0; i<n_models; i++){

			temp[i] = (double)(std::rand() % 100 + 1);					// generate number between 1 and 100
			sum += temp[i];
		}

		for (int i=0; i<n_models; i++)
			Rni(n,i) = temp[i]/sum;

		delete [] temp; temp = NULL;

	}
	initialize_helper();									// Initialize other parameters.
}

void PPCA_Mixture_EM::initialize_kmeans(int f_max_iter){

	kMeans classifier(n_models, f_max_iter, data);			// cluster data

	classifier.initCentersRandom();

	std::vector <int> groups = classifier.findSolution();

	mean = classifier.getCenters();							// Initialize mean

	int total = n_obs;									// Initialize model priors.

	for (int n=0; n<n_obs; n++){

		mixfrac[groups[n]] += 1.0;
		total ++;
	}
	/**
	for (int i=0; i<n_models; i++){						// make sure that no priors are 0
		if (mixfrac[i]==0){
			mixfrac[i] = 1;
			total++;
		}
	}
	**/
	for (int i=0; i<n_models; i++)
		mixfrac[i] /= (double)total;

	int * total_per_model = new int[n_models]();

	for (int n=0; n<n_obs; n++){						// Initialize data covariance matrix with data points assigned to each model by kmean.

		vec temp = data.col(n) - mean.col(groups[n]);

		S_mat_vector[groups[n]] += temp*temp.t();

		total_per_model[groups[n]] ++;

	}

	for (int i=0; i<n_models; i++){								// TODO: if we get a bad start, re-initialize.

		std::cout << "Initial data points assigned to model " << i << ": "
				<<total_per_model[i] <<"\n";					// Output to screen the number of data points assigned to each model.
	}
	std::cout << "\n";

	for (int i=0; i<n_models; i++)
		S_mat_vector[i] /= (double)total_per_model[i];

	delete[] total_per_model; total_per_model = NULL;


	for (int i=0; i<n_models; i++){

		mat Si_eigvec;
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


		mat tempU; mat tempV; vec temps(n_components); vec temps_inv(n_components);

		svd(tempU, temps, tempV, tempM);

		for (int n=0; n<n_components; n++){
			if (std::abs(temps(n)) < SMALL_NUMBER){
				temps_inv(n) = 0;
			}
			else
				temps_inv(n) = 1.0 / temps(n);
		}

		Minv_mat_vector[i] = tempV * diagmat(temps_inv) * tempU.t();
	}

	write_to_file_mean(0);									// Write to file itial mean.
	write_to_file_params(0);

	has_init = 1;
}


void PPCA_Mixture_EM::initialize_helper(void){

	// The initialize random functions above only assign values to Rni. From this, the rest of the parameters are updated.

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
		mat Si_eigvec;
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

	write_to_file_mean(0);									// Write to file initial mean.
	write_to_file_params(0);
	// All parameters initialized, ready to start iterating.
	has_init = 1;

}

void PPCA_Mixture_EM::update_Rni_all( void){

	double * total = new double[n_obs]();
	double * max = new double[n_obs]();		// Store largest log_Ptn_i for each observation. This will be subtracted for each log_Ptn_i to avoid overflow.

	mat log_Ptn_i(n_obs,n_models,fill::zeros);

	for (int i=0; i<n_models; i++){

		mat C = eye<mat>(n_var,n_var) * noise_var[i] +  W_mat_vector[i] * W_mat_vector[i].t();

		//mat Cinv = (eye<mat>(n_var,n_var) - W_mat_vector[i] * Minv_mat_vector[i] * W_mat_vector[i].t()) / noise_var[i];
		mat Cinv;

		mat U_C; mat V_C; vec s_C(n_var); vec s_C_inv(n_var);

		svd(U_C,s_C,V_C,C);

		double log_det_C(0);
		for (int k=0; k<n_var; k++){
			if (std::abs(s_C(k)) < SMALL_NUMBER){
				//log_det_C += std::log(SMALL_NUMBER);
				s_C_inv(k) = 0;
			}
			else{
				if (s_C(k) < 0)
					log_det_C += std::log(- s_C(k));
				else
					log_det_C += std::log(s_C(k));
				s_C_inv(k) = 1.0/s_C(k);
			}
		}

		Cinv = V_C * diagmat(s_C_inv) * U_C.t();

		for (int n=0; n<n_obs; n++)
			log_Ptn_i(n,i) = calc_log_Ptn_i(n,i,Cinv, log_det_C);
	}

		for (int n=0; n<n_obs; n++)
			max[n] = log_Ptn_i.row(n).max();

		for (int n=0; n<n_obs; n++){								// subtract largest log_Ptn_i from each model for each observation

			for (int i=0; i<n_models; i++){							// This is all for avoiding large numbers. Better numerical stability.

				log_Ptn_i(n,i) -= (max[n]) ;

				total[n] += std::exp(log_Ptn_i(n,i)) * mixfrac[i];

				#ifdef DEBUG
					if (isinf(total[n]) )
						std::cout << "temp_total for obs " << n << " is inf\n";
				#endif
			}
		}

		log_likelihood = 0;											// Calculate log likelihood.

		for (int n=0; n<n_obs; n++)
			log_likelihood += std::log(total[n]);

	for (int n=0; n<n_obs; n++){									// Normalize.
		for(int i=0; i<n_models; i++){

			Rni(n,i) = std::exp(log_Ptn_i(n,i)) / total[n] * mixfrac[i] ;

			if (Rni(n,i) <= SMALL_NUMBER )	// Set lower limit to Rni to avoid mix_frac from being zero in the case that all obs for this n are 0.
				Rni(n,i) = SMALL_NUMBER;

			#ifdef DEBUG
				if (n ==0){
					std::cout << "argument to Rni: " << log_Ptn_i(n,i) + std::log(mixfrac[i]) - std::log(total[n])<< "\n";
					std::cout << "std::exp(log_Pt_0_" << i << "): " << log_Ptn_i(n,i) << "\n";
					if (Rni(n,i) == SMALL_NUMBER)
						std::cout << "Rn(0," << i  << ") was floored to SMALL_NUMBER: " << SMALL_NUMBER <<"\n\n";
					else
						std::cout << "Rn(0," << i  << "): " <<  Rni(n,i) << "\n\n";
				}
				if (std::isinf(Rni(n,i)))
					std::cout << "INF detected Rni_i\n";
				if (std::isnan(Rni(n,i))){
					std::cout << "NAN detected Rn_i_2\n";
					std::cout << "log_Ptn_i(" << n <<","<<i<<"): "<<log_Ptn_i(n,i) << "\n";


					std::cout << "Total[" << n << "]: " <<total[n] << "\n";
					std::cout << "mixfrac[" << i << "]: " << mixfrac[i] << "\n";
				}
			#endif
		}
	}

	delete [] total; total = NULL;
	delete [] max; max = NULL;
}

double PPCA_Mixture_EM::calc_log_Ptn_i(int f_n, int f_i, mat f_Cinv, double f_log_det_C){
	colvec temp = data.col(f_n) - mean.col(f_i);

	double arg = as_scalar(temp.t()* f_Cinv * temp) / 2.0;

	double log_Ptn_i = -0.5 * f_log_det_C -arg;

	return log_Ptn_i;
}

double PPCA_Mixture_EM::calc_log_Ptn_i(mat &f_samples, int f_n, int f_i, mat f_Cinv, double f_log_det_C){
	colvec temp = f_samples.col(f_n) - mean.col(f_i);

	double arg = as_scalar(temp.t()* f_Cinv * temp) / 2.0;

	double log_Ptn_i = -0.5 * f_log_det_C -arg;

	return log_Ptn_i;
}

double PPCA_Mixture_EM::calc_log_Ptn_i(mat &f_samples, mat &f_mean, int f_n, int f_i, mat f_Cinv, double f_log_det_C){
	colvec temp = f_samples.col(f_n) - f_mean.col(f_i);

	double arg = as_scalar(temp.t()* f_Cinv * temp) / 2.0;

	double log_Ptn_i = -0.5 * f_log_det_C -arg;

	return log_Ptn_i;
}


void PPCA_Mixture_EM::optimize(int f_max_iter, int write_freq_mean, int write_freq_params){

	double MIN_NOIS_VAR = SMALL_NUMBER;
	std::vector<double> init_noise_var = noise_var;	// Store initial noise_var to avoid small values.

	if (!has_init)									// Throw exception if calling this function before initializing parameters.
		throw no_init_exception();

	n_iter = 0;										// reset

	while (n_iter < f_max_iter){

		update_Rni_all();

		double total_mixfrac(0);

		for (int i=0; i<n_models; i++){
			mixfrac[i] = 0.0;								// Update model prior probability.

			for (int n=0; n<n_obs; n++)
				mixfrac[i] += Rni(n,i);

			mixfrac[i] /= (double)n_obs;

			total_mixfrac += mixfrac[i];

			#ifdef DEBUG
				std::cout << "mixfrac[" << i << "]: " << mixfrac[i] << "\n";
   	   	   	#endif
		}

		#ifdef DEBUG
			std::cout << "total_mixfrac: " << total_mixfrac;
			std::cout << "\n";
	   	#endif

		for (int i=0; i<n_models; i++){

			colvec temp_mean(n_var,fill::zeros);			// Update mean parameter.
			double temp_den(0);

			for (int n=0; n<n_obs; n++){

				temp_mean += Rni(n,i) * data.col(n);
				temp_den += Rni(n,i);
			}

			mean.col(i) = temp_mean / temp_den;

			#ifdef DEBUG
				if (mean.col(i).has_nan())
					std::cout << "NAN detected: mean for model " << i << "\n";
			#endif
		}

		for (int i=0; i<n_models; i++){						// Update W mat, noise_var, S mat, and Minv mat last

			mat temp_Si(n_var,n_var, fill::zeros);									// Update S matrices.

			for (int n=0; n<n_obs; n++){

				colvec tempvec = data.col(n) - mean.col(i);
				temp_Si += Rni(n,i)*tempvec*tempvec.t();
			}

			S_mat_vector[i] = temp_Si/(mixfrac[i] * n_obs);

			mat tempSW_mat(n_var, n_components, fill::zeros);			//Calculate S*W.
			/**
			 * The below lines are meant to calculate S*W efficiently, but we need to divide by
			 * mixfrac[i]*n_obs.
			 * **/
			//for (int n=0; n<n_obs; n++){
			//
			//	colvec tempvec = data.col(n) - mean(i);
			//	tempSW_mat += Rni(n,i) * tempvec * (tempvec.t() * W_mat_vector[i]);
			//}
			//tempSW_mat /= (mixfrac[i]*n_obs);
			tempSW_mat =  S_mat_vector[i] *W_mat_vector[i];

			mat tempInv_U; mat tempInv_V; vec tempinvs(n_components); vec tempinvs_inv(n_components);

			svd(tempInv_U, tempinvs, tempInv_V, noise_var[i] * eye(n_components, n_components) + Minv_mat_vector[i] * W_mat_vector[i].t() * tempSW_mat);

			for (int n=0; n<n_components; n++){
				if (std::abs(tempinvs(n)) < SMALL_NUMBER){
					tempinvs_inv(n) = 0;

					#ifdef DEBUG
						std::cout << "temps n: " <<n << " i: " << i <<"\n";
					#endif
				}
				else
					tempinvs_inv(n) = 1.0 / tempinvs(n);
			}

			mat tempInv_mat = tempInv_V * diagmat(tempinvs_inv) * tempInv_U.t();

			//mat tempInv_mat = inv(noise_var[i] * eye(n_components, n_components) + Minv_mat_vector[i] * W_mat_vector[i].t() * tempSW_mat);

			W_mat_vector[i] = tempSW_mat * tempInv_mat;					// Update W mat.

			noise_var[i] = 1.0/((double)n_var) * trace(S_mat_vector[i] - tempSW_mat * Minv_mat_vector[i] * W_mat_vector[i].t());

			if (noise_var[i]<MIN_NOIS_VAR)							// Check for small noise variance values and replace with initial value if too small.
				noise_var[i] = init_noise_var[i];

			//mat tempM = eye<mat>(n_components, n_components);


			mat tempU; mat tempV; vec temps(n_components); vec temps_inv(n_components);

			svd(tempU, temps, tempV, noise_var[i] * eye<mat>(n_components, n_components) + W_mat_vector[i].t() * W_mat_vector[i]);

			for (int n=0; n<n_components; n++){
				if (std::abs(temps(n)) < SMALL_NUMBER){
					temps_inv(n) = 0;

					#ifdef DEBUG
						std::cout << "temps n: " <<n << " i: " << i <<"\n";
					#endif
				}
				else
					temps_inv(n) = 1.0 / temps(n);
			}


			Minv_mat_vector[i] = tempV * diagmat(temps_inv) * tempU.t();

			//Minv_mat_vector[i] = inv(noise_var[i] * eye<mat>(n_components, n_components) + W_mat_vector[i].t() * W_mat_vector[i]); // Update Minv matrices.
		}

		std::cout << "Iteration " << n_iter + 1 << " completed.\n";

		this->n_iter++;

		if ((n_iter % write_freq_mean)==0)
			write_to_file_mean(n_iter);									// Write to file mean.
		if ((n_iter % write_freq_params)==0)						// Write to file params.
			write_to_file_params(n_iter);

	}
}

PPCA_Mixture_SAG::PPCA_Mixture_SAG(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, int f_n_mini_batch):

	PPCA_Mixture_EM(f_data, f_dim, f_n_components, f_n_models),

	mixfrac_softmax_coeff(f_n_models,0),  n_mini_batch(f_n_mini_batch), lip_const(0), base_rate(0), grad_mean(f_dim[0]*f_dim[1],f_n_models,fill::zeros), grad_noise_var(f_n_models,0),

	grad_mixfrac_softmax_coeff(f_n_models,0), grad_W_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros)),

	grad_mean_prev(f_dim[0]*f_dim[1],f_n_models,fill::zeros), grad_noise_var_prev(f_n_models,0), grad_mixfrac_softmax_coeff_prev(f_n_models,0),

	grad_W_mat_vector_prev(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros))
{

}


PPCA_Mixture_SAG::PPCA_Mixture_SAG(mat f_data, std::vector<int> f_dim, int f_n_components, int f_n_models, char* prev_run_filename, int f_n_mini_batch):

	PPCA_Mixture_EM(f_data, f_dim, f_n_components, f_n_models, prev_run_filename),

	mixfrac_softmax_coeff(f_n_models,0),  n_mini_batch(f_n_mini_batch), lip_const(0), base_rate(0), grad_mean(f_dim[0]*f_dim[1],f_n_models,fill::zeros), grad_noise_var(f_n_models,0),

	grad_mixfrac_softmax_coeff(f_n_models,0), grad_W_mat_vector(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros)),

	grad_mean_prev(f_dim[0]*f_dim[1],f_n_models,fill::zeros), grad_noise_var_prev(f_n_models,0), grad_mixfrac_softmax_coeff_prev(f_n_models,0),

	grad_W_mat_vector_prev(f_n_models,mat(f_dim[0]*f_dim[1],f_n_components,fill::zeros))
{

}

void PPCA_Mixture_SAG::initialize_random_SAG(void){
	initialize_random();
	base_rate = 1.0/16.0;
	lip_const = 1.0;
	mixfrac_softmax_coeff = mixfrac_to_softmax(this->mixfrac);

}

void PPCA_Mixture_SAG::initialize_kmeans_SAG(int f_max_iter){
	initialize_kmeans(f_max_iter);
	base_rate = 1.0/16;
	lip_const = 1.0;
	mixfrac_softmax_coeff = mixfrac_to_softmax(this->mixfrac);
}

std::vector<double> PPCA_Mixture_SAG::mixfrac_to_softmax(std::vector<double> f_mix_frac_coeff){

	std::vector<double> softmax_values(f_mix_frac_coeff.size(),0);

	for (int i=0; i<f_mix_frac_coeff.size(); i++){
		softmax_values[i] = std::log(f_mix_frac_coeff[i]);
	}
	return softmax_values;
}

std::vector<double> PPCA_Mixture_SAG::softmax_to_mixfrac(std::vector<double> f_soft_max_coeff){

	double norm(0);
	std::vector<double> mix_frac_values(f_soft_max_coeff.size(),0);

	for (int i=0; i<f_soft_max_coeff.size(); i++){
		mix_frac_values[i] = std::exp(f_soft_max_coeff[i]);
		norm += mix_frac_values[i];
	}
	for (int i=0; i<f_soft_max_coeff.size(); i++){
		mix_frac_values[i] /= norm;
	}

	return mix_frac_values;
}

mat PPCA_Mixture_SAG::get_mini_batch_SAG(int f_n_mini_batch, double max_trans){

	mat outmat(data_dim[0]*data_dim[1], f_n_mini_batch);

	std::srand(std::time(NULL));							// generate seed.

	for (int i=0; i<f_n_mini_batch; i++){

		int indx = std::rand() % n_obs;
		double * temp = new double[data_dim[0]*data_dim[1]];

		double d, theta, delta_x, delta_y;

		d = (double)rand() / RAND_MAX;
		theta = d * 360.0;

		d = (double)rand() / RAND_MAX;
		delta_x = -max_trans + d * (max_trans*2);

		d = (double)rand() / RAND_MAX;
		delta_y = -max_trans + d * (max_trans*2);

		vec tcol = data.col(indx);
		transform_Img_single(tcol.memptr(), temp, data_dim[0], data_dim[1], delta_x, delta_y, theta, (double)std::min(data_dim[0],data_dim[1]));

		for (int k=0; k<data_dim[0]*data_dim[1]; k++){
			outmat(k,indx) = temp[k];
		}
		delete [] temp; temp = NULL;
	}

	return outmat;
}

double PPCA_Mixture_SAG::calc_base_rate_SAG(void){

	return std::max(1.0/16, std::pow(2,1.0-(double)n_iter/150.0));

}

std::vector<double> PPCA_Mixture_SAG::calc_Cinv_log_det_C(std::vector<mat> &f_W_mat_vector, std::vector<double> f_noise_var, std::vector<mat> &f_out_mat_vector){

	std::vector<double> log_det_C_vector(n_models,0);

	for (int i=0; i<n_models; i++){

		mat C = eye<mat>(n_var,n_var) * f_noise_var[i] +  f_W_mat_vector[i] * f_W_mat_vector[i].t();

		mat U_C; mat V_C; vec s_C(n_var); vec s_C_inv(n_var,fill::zeros);

		svd(U_C,s_C,V_C,C);

		for (int k=0; k<n_var; k++){
			if (std::abs(s_C(k)) < SMALL_NUMBER){
				s_C_inv(k) = 0;
			}
			else{
				if (s_C(k) < 0)
					log_det_C_vector[i] += std::log(- s_C(k));
				else
					log_det_C_vector[i] += std::log(s_C(k));
				s_C_inv(k) = 1.0/s_C(k);

			}
		}

		f_out_mat_vector[i] = V_C * diagmat(s_C_inv) * U_C.t();
	}
	return log_det_C_vector;
}

void PPCA_Mixture_SAG::calc_log_Ptn_i_mat(std::vector<mat> &f_Cinv_vector, std::vector<double> f_log_det_C_vector, mat &f_mean, mat &f_samples, int f_n_samples, mat &f_out_mat){

	for (int i=0; i<n_models; i++){

		for (int n=0; n<f_n_samples; n++)
			f_out_mat(n,i) = calc_log_Ptn_i(f_samples, f_mean, n, i, f_Cinv_vector[i], f_log_det_C_vector[i]);
	}

	double * max = new double[f_n_samples]();

	for (int n=0; n<f_n_samples; n++){
		max[n] = f_out_mat.row(n).max();

		for (int i=0; i<n_models; i++)							// This is all for avoiding large numbers. Better numerical stability.
			f_out_mat(n,i) -= (max[n]) ;
	}
}

void PPCA_Mixture_SAG::calc_log_Ptn_i_mat_no_norm(std::vector<mat> &f_Cinv_vector, std::vector<double> f_log_det_C_vector, mat &f_mean, mat &f_samples, int f_n_samples, mat &f_out_mat){
	for (int i=0; i<n_models; i++){

		for (int n=0; n<f_n_samples; n++)
			f_out_mat(n,i) = calc_log_Ptn_i(f_samples, f_mean, n, i, f_Cinv_vector[i], f_log_det_C_vector[i]) - 0.5 * n_var * std::log(2*3.14);
	}
}

double PPCA_Mixture_SAG::estimate_log_likelihood(mat &f_log_Ptn_i_mat, std::vector<double> f_mixfrac, mat &f_samples, int f_n_samples){

	double * total = new double[f_n_samples]();

	for (int n=0; n<f_n_samples; n++){

		for (int i=0; i<n_models; i++){
			total[n] += std::exp(f_log_Ptn_i_mat(n,i)) * f_mixfrac[i];
		}

	}

	double out_likelihood(0);

	for (int n=0; n<f_n_samples; n++)
		out_likelihood += std::log(total[n]);

	return out_likelihood / f_n_samples;
}

double PPCA_Mixture_SAG::estimate_Lip_SAG(int iter, double f_k_1, double f_k, double f_grad_k){

	if (n_iter%iter){
		return (double)n_iter / std::pow(2,1.0/150.0);
	}

	else {
		//double L(1);
		//while (f_k_1 - f_k <= f_grad_k)
	}
	return 0.0;				// TODO implements this...
}

std::vector<double> PPCA_Mixture_SAG::calc_grad_SAG(std::vector<mat> &f_Cinv_vector, mat &f_log_Ptn_i_mat,std::vector<mat> &f_W_mat_vector, mat &f_mean, std::vector<double> f_mix_frac,
		std::vector<double> f_mix_frac_softmax_coeff, std::vector<double> f_noise_var, mat &f_samples, int f_n_samples ){

	std::vector<double> norm(f_n_samples,0);
	double  norm_mix_frac_softmax_coeff(0);

	mat temp_grad_mean(n_var,n_models,fill::zeros);

	std::vector<double> temp_grad_noise_var(n_models,0);

	std::vector<double> temp_grad_mixfrac(n_models,0);

	std::vector<double> temp_grad_mix_frac_softmax_coeff(n_models,0);

	std::vector<mat> temp_grad_C_mat_vector(n_models,mat(n_var,n_var,fill::zeros));

	std::vector<mat> temp_grad_W_mat_vector(n_models,mat(data_dim[0]*data_dim[1],n_components,fill::zeros));


	for (int i=0; i<n_models; i++){																	// Pre-calculate the normalization constant for each particle.

		norm_mix_frac_softmax_coeff += std::exp(f_mix_frac_softmax_coeff[i]);

		for (int n=0; n<f_n_samples; n++){
			norm[n] += f_mix_frac[i] * std::exp(f_log_Ptn_i_mat(n,i));
		}
	}

	for (int i=0; i<n_models; i++){
		for (int n=0; n<f_n_samples; n++){

			double numer = f_mix_frac[i] * std::exp(f_log_Ptn_i_mat(n,i));

			vec temp_vec = f_Cinv_vector[i] * (f_samples.col(n)-f_mean.col(i));						// Pre-calculate Cinv * (tn-ui) used for calculating gradients.

			temp_grad_mean.col(i) += (numer / norm[n] * temp_vec);

			temp_grad_mixfrac[i] += (numer / f_mix_frac[i] / norm[n]) ;

			for (int j=0; j<n_models; j++){
				if (i==j)
					temp_grad_mix_frac_softmax_coeff[i]  += (std::exp(f_log_Ptn_i_mat(n,j)) / norm[n] *
							(std::exp(f_mix_frac_softmax_coeff[i])/norm_mix_frac_softmax_coeff
								- std::pow(std::exp(f_mix_frac_softmax_coeff[i]),2)/std::pow(norm_mix_frac_softmax_coeff,2))) ;
				else
					temp_grad_mix_frac_softmax_coeff[i]  += (std::exp(f_log_Ptn_i_mat(n,j)) / norm[n] * (- std::exp(f_mix_frac_softmax_coeff[i]) *
							std::exp(f_mix_frac_softmax_coeff[j]) / std::pow(norm_mix_frac_softmax_coeff,2)));
			}

			temp_grad_noise_var[i]  += (numer /norm[n] * 0.5 * (-trace(f_Cinv_vector[i]) + as_scalar(temp_vec.t() * temp_vec))); //This calculation is the derivative with respect to sig^2, and does not constrain noise_var to be positive.

			temp_grad_C_mat_vector[i] += (numer / norm[n] * 0.5 * (-f_Cinv_vector[i] + temp_vec * temp_vec.t()));

		}

		temp_grad_mixfrac[i] /= (double)f_n_samples;
		temp_grad_mix_frac_softmax_coeff[i] /= (double)f_n_samples;
		temp_grad_noise_var[i] /= (double)f_n_samples;
		temp_grad_mean.col(i) /= (double)f_n_samples;
	}

	/**
	for (int i=0; i<n_models; i++){								// Use chain rule to calculate gradients with respect to softmax coefficient.

		for (int j=0; j<n_models; j++){

			if (i==j)
				temp_grad_mix_frac_softmax_coeff[i] += temp_grad_mixfrac[i] * (std::exp(f_mix_frac_softmax_coeff[i])/norm_mix_frac_softmax_coeff[i] - std::pow(std::exp(f_mix_frac_softmax_coeff[i]),2)/std::pow(norm_mix_frac_softmax_coeff[i],2));
			else
				temp_grad_mix_frac_softmax_coeff[i] += temp_grad_mixfrac[i] * (- std::exp(f_mix_frac_softmax_coeff[i]) * std::exp(f_mix_frac_softmax_coeff[j]) / std::pow(norm_mix_frac_softmax_coeff[i],2));
		}

	}
	**/
	/**
	for (int i=0; i<n_models; i++){								// Use chain rule to calculate gradients with respect to W mat.
		for (int j=0; j<data_dim[0]*data_dim[1]; j++){
			for (int k=0; k<n_components; k++){

				mat temp_mat(n_var,n_var,fill::zeros);

				temp_mat.col(j) += f_W_mat_vector[i].col(k);

				temp_mat.row(j) += f_W_mat_vector[i].col(k).t();

				temp_grad_W_mat_vector[i](j,k) = trace(temp_grad_C_mat_vector[i].t() * temp_mat);
			}
		}
	}
	// Copy all gradient data into a vector. This is inefficient, but to make code a bit clearer. Should try to avoid this step in the future.
	std::vector<double> out_vector( temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size() + temp_grad_mix_frac_softmax_coeff.size() + n_models*temp_grad_W_mat_vector[0].n_rows*temp_grad_W_mat_vector[0].n_cols,0);

	for (int k=0; k<temp_grad_mean.n_rows*temp_grad_mean.n_cols; k++)
		out_vector[k] = temp_grad_mean(k);

	for (int k=0; k<temp_grad_noise_var.size(); k++)
		out_vector[k + temp_grad_mean.n_rows*temp_grad_mean.n_cols] = temp_grad_noise_var[k];

	for (int k=0; k<temp_grad_mix_frac_softmax_coeff.size(); k++)
		out_vector[k + temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size()] = temp_grad_mix_frac_softmax_coeff[k];

	for (int k=0; k<n_models ; k++){
		int num_elements = temp_grad_W_mat_vector[0].n_rows * temp_grad_W_mat_vector[0].n_cols;
		for (int j=0; j<num_elements; j++)
			out_vector[k*num_elements + j + temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size() + temp_grad_mix_frac_softmax_coeff.size()] = temp_grad_W_mat_vector[k](j);
	}
	**/

	std::vector<double> out_vector(n_models);
	for (int k=0; k<n_models; k++)
		out_vector[k] = temp_grad_mixfrac[k];
	return out_vector;

}

std::vector<double> PPCA_Mixture_SAG::calc_grad_finite_dif_SAG(std::vector<mat> f_W_mat_vector, mat f_mean, std::vector<double> f_mix_frac,
		std::vector<double> f_mix_frac_softmax_coeff, std::vector<double> f_noise_var, mat &f_samples, int f_n_samples ){

	double esp = 10e-5;

	std::vector<double> norm(f_n_samples,0);

	std::vector<double> norm_mix_frac_softmax_coeff(n_models,0);

	std::vector<mat> temp_Cinv_mat_vector(n_models, mat(n_var,n_var,fill::zeros));

	mat temp_grad_mean(n_var,n_models,fill::zeros);

	std::vector<double> temp_grad_noise_var(n_models,0);

	std::vector<double> temp_grad_mixfrac(n_models,0);

	std::vector<double> temp_grad_mix_frac_softmax_coeff(n_models,0);

	std::vector<mat> temp_grad_W_mat_vector(n_models,mat(data_dim[0]*data_dim[1],n_components,fill::zeros));

	std::vector<double> temp_log_det_C = calc_Cinv_log_det_C(f_W_mat_vector, f_noise_var, temp_Cinv_mat_vector);

	mat temp_log_Ptn_i_mat(f_n_samples,n_models);

	// Calculate gradients for means.
	for (int i=0; i<n_models; i++){
		for (int j=0; j<n_var; j++){

			f_mean.col(i)(j) += esp;																		// Perturb a single element.

			calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);

			double f1 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

			f_mean.col(i)(j) -=(2*esp);

			calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);

			double f2 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

			f_mean.col(i)(j) += esp;

			temp_grad_mean.col(i)(j) = (f1-f2) / (2*esp);
		}
	}
	//Calculate gradiets for mixfrac.
	for (int i=0; i<n_models; i++){

			f_mix_frac[i] += esp;																		// Perturb a single element.

			calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);
			double f1 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

			f_mix_frac[i] -= (2*esp);

			calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);
			double f2 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

			f_mix_frac[i] += esp;

			temp_grad_mixfrac[i] = (f1-f2) / (2*esp);


	}

	// Calculate gradients for mix_frac_soft_max.
	for (int i=0; i<n_models; i++){

		f_mix_frac_softmax_coeff[i] += esp;
		f_mix_frac  = softmax_to_mixfrac(f_mix_frac_softmax_coeff);

		calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);
		double f1 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

		f_mix_frac_softmax_coeff[i] -= (2*esp);
		f_mix_frac = softmax_to_mixfrac(f_mix_frac_softmax_coeff);

		calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);
		double f2 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

		f_mix_frac_softmax_coeff[i] += esp;
		f_mix_frac = softmax_to_mixfrac(f_mix_frac_softmax_coeff);

		temp_grad_mix_frac_softmax_coeff[i] = (f1-f2)/ (2*esp);

	}

	// Calculate gradients for noise var.
	for (int i=0; i<n_models; i++){

		f_noise_var[i] += esp;

		temp_log_det_C = calc_Cinv_log_det_C(f_W_mat_vector, f_noise_var, temp_Cinv_mat_vector);
		calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);

		double f1 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

		f_noise_var[i] -= (2*esp);

		temp_log_det_C = calc_Cinv_log_det_C(f_W_mat_vector, f_noise_var, temp_Cinv_mat_vector);
		calc_log_Ptn_i_mat_no_norm(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);

		double f2 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

		f_noise_var[i] += esp;

		temp_grad_noise_var[i] = (f1-f2)/(2*esp);

	}
	/**
	// Calculate gradients for W mat.
	for (int i=0; i<n_models; i++){
		for (int j=0; j<n_var; j++){
			for (int k=0; k<n_models; k++){

				f_W_mat_vector[i](j,k) += esp;

				temp_log_det_C = calc_Cinv_log_det_C(f_W_mat_vector, f_noise_var, temp_Cinv_mat_vector);
				calc_log_Ptn_i_mat(temp_Cinv_mat_vector, temp_log_det_C, f_mean, f_samples, f_n_samples, temp_log_Ptn_i_mat);
				double f1 = estimate_log_likelihood(temp_log_Ptn_i_mat, f_mix_frac, f_samples, f_n_samples);

				f_W_mat_vector[i](j,k) -= esp;

				temp_grad_W_mat_vector[i](j,k) = (f1-f0)/esp;

			}
		}

	}


	// Copy all gradient data into a vector. This is inefficient, but to make code a bit clearer. Should try to avoid this step in the future.
	std::vector<double> out_vector( temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size() + temp_grad_mix_frac_softmax_coeff.size() + n_models*temp_grad_W_mat_vector[0].n_rows*temp_grad_W_mat_vector[0].n_cols,0);

	for (int k=0; k<temp_grad_mean.n_rows*temp_grad_mean.n_cols; k++)
		out_vector[k] = temp_grad_mean(k);

	for (int k=0; k<temp_grad_noise_var.size(); k++)
		out_vector[k + temp_grad_mean.n_rows*temp_grad_mean.n_cols] = temp_grad_noise_var[k];

	for (int k=0; k<temp_grad_mix_frac_softmax_coeff.size(); k++)
		out_vector[k + temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size()] = temp_grad_mix_frac_softmax_coeff[k];

	for (int k=0; k<n_models ; k++){
		int num_elements = temp_grad_W_mat_vector[0].n_rows * temp_grad_W_mat_vector[0].n_cols;
		for (int j=0; j<num_elements; j++)
			out_vector[k*num_elements + j + temp_grad_mean.n_rows*temp_grad_mean.n_cols + temp_grad_noise_var.size() + temp_grad_mix_frac_softmax_coeff.size()] = temp_grad_W_mat_vector[k](j);
	}
	**/
	std::vector<double> out_vector(n_models);
	for (int k=0; k<n_models; k++)
		out_vector[k] = temp_grad_mixfrac[k];
	return out_vector;
}

