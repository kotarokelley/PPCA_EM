/*
 * Exceptions.h
 *
 *  Created on: Sep 15, 2016
 *      Author: kotarokelley
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

/***---Exceptions---***/
class matrix_compatibility_exception: public std::exception {
	public:
		void err_msg(void);
};

class no_init_exception: public std::exception {
	public:
		void err_msg(void);
};

#endif /* EXCEPTIONS_H_ */
