/*
 * Exceptions.cpp
 *
 *  Created on: Dec 10, 2016
 *      Author: kotarokelley
 */

#include "Exceptions.h"


void matrix_compatibility_exception::err_msr(void){
	printf("\nMatrix dimensions did not match. Nothing was modified\n\n");
} //mat_err;


void no_init_exception::err_msg(void) {
	printf("\nA call to an initialization function must be made before EM optimization\n\n");
} //init_err;
