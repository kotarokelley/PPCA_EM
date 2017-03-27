/*
 * Exceptions.h
 * Author: kotarokelley
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***---Exceptions---***/
class matrix_compatibility_exception: public std::exception {
	public:
		void err_msg(void);
};

class no_init_exception: public std::exception {
	public:
		void err_msg(void);
};

class image_operations_exception: public std::exception {
	public:
		image_operations_exception(int f_function);
		int function;
		void err_msg(void);
};

#endif /* EXCEPTIONS_H_ */
