
#include "class2D.h"

using namespace arma;

/**---Class Constructor---***/
class2D::class2D(float* f_data, int* f_dim, int f_n_classes):

	data(f_data,f_dim[0]*f_dim[1],f_dim[2]), n_classes(f_n_classes), n_components(10) {
	data_dim[0] = f_dim[0]*f_dim[1];
	data_dim[1] = f_dim[2];
}

class2D::class2D(float* f_data, int* f_dim, int f_n_classes, int f_n_components):

	data(f_data,f_dim[0]*f_dim[1],f_dim[2]), n_classes(f_n_classes), n_components(f_n_components) {
	data_dim[0] = f_dim[0]*f_dim[1];
	data_dim[1] = f_dim[2];
}

/**---Class Methods---***/
void class2D::classify_PPCA_EM(double f_trans, double f_trans_rate, double f_ang, double f_ang_rate, int f_max_iter){

	//PPCA_Mixture_EM* pca = new PPCA_Mixture_EM(this->data, this->data_dim);



	//std::cout << pca->get_params();
	 //Clean up.
	//delete pca; pca = NULL;

}
void class2D::tofile(char * f_name){

}



