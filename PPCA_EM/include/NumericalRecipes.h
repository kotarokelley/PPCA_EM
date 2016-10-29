/*
 * NumericalRecipes.h
 *
 *  Created on: Sep 13, 2016
 *      Author: kotarokelley
 */

#ifndef NUMERICALRECIPES_H
#define NUMERICALRECIPES_H

struct Mat;

/***---Functions---***/
int matrixMult( Mat *, Mat *, Mat *);
/*	Arguments:
 * 		Pointers to each in A*B = C.
 * 		The sizes of each matrix must be correct and will be checed only through the rows and cols fields in each Mat struct.
 * 	Returns:
 * 		If success, returns 1.
 * 		If fail due to size incompatibility, returns 0.
 * 		If fial due to some other reason, returns -1.
 * 	Exceptions:
 *		out_of_range
 */
void printMatrix( Mat *);

#endif /* NUMERICALRECIPES_H */
