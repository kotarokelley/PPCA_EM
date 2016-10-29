/*
 * DataStructures.h
 *
 *  Created on: Sep 15, 2016
 *      Author: kotarokelley
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

/***---Data Structures---***/
struct Mat {
	/**	Matrix data structure.
	 * 	A 2-D matrix that is stored as a 1-D array.
	 *
	 * 	Class Methods:
	 *
	 * 	Class Members:
	 * 		rows: int
	 * 			Number of rows.
	 * 		cols: int
	 * 			Number of columns.
	 */
	int rows;
	int cols;
	float * dat;

	Mat(int f_rows, int f_cols) : rows(f_rows), cols(f_cols) {
		dat = new float[rows*cols];
	}

	Mat(int f_rows, int f_cols, float *f_dat) : rows(f_rows), cols(f_cols) {
		dat = new float[rows*cols];
		for (int i=0; i<rows*cols; i++){
			dat[i] = f_dat[i];
		}
	}

	~Mat(){
		delete [] dat;
		dat = NULL;
	}
};

#endif /* DATASTRUCTURES_H_ */
