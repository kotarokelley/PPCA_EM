//============================================================================
// Name        : NumericalRecipes.cpp
// Author      : Kotaro Kelley
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <exception>
#include <limits.h>
#include <math.h>
#include "NumericalRecipes.h"
#include "Exceptions.h"
#include "DataStructures.h"

using namespace std;

/***---Functions---***/
int matrixMult(  Mat * A, Mat * B, Mat * C){

	try{
		if ( A->cols != B->rows || A->rows != C->rows || B->cols !=C->cols)
			throw mat_err;
		else {
			for (int i=0; i<A->rows; i++)
				for (int j=0; j<B->cols; j++)
					for (int k=0; k<A->cols; k++)
						C->dat[i*C->cols+j] = C->dat[i*C->cols+j] + A->dat[i*A->cols+k]*B->dat[k*B->cols+j];
		return 1;
		}
	}
	catch(matrix_compatibility_exception& e){
		//e.err_msg();
		return 0;
	}
	catch(std::exception& e){
		//printf("Something went wrong with matrix multiplication.\n");
		return -1;
	}

}





void printMatrix(Mat * mat){
	for (int i=0; i<mat->rows; i++){
		for (int j=0; j<mat->cols-1; j++){
			printf("%5e ", mat->dat[mat->cols*i + j]);
		}
		printf("%5e \n",mat->dat[mat->cols*i + mat->cols-1]);
	}
}

