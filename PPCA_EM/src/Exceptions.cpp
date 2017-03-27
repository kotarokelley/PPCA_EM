/*
 * Exceptions.cpp
 *
 *      Author: kotarokelley
 */

#include "Exceptions.h"


void matrix_compatibility_exception::err_msg(void){
	std::cout << "\nMatrix dimensions did not match. Nothing was modified\n\n";
} //mat_err;


void no_init_exception::err_msg(void) {
	std::cout << "\nA call to an initialization function must be made before EM optimization\n\n";
} //init_err;


image_operations_exception::image_operations_exception(int f_function): function(f_function){}

void image_operations_exception::err_msg(void){

	std::cout << "\nA call to ImageOperations function: ";

	switch(this->function){
		case 0:
			std::cout << "resample_Img\n";
			break;
		case 1:
			std::cout << "rebin_Img\n";
			break;
		default:
			std::cout << "unknown function\n";
			break;
	}
}
