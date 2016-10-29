/*
 * test_NumericalRecipes.cpp
 *
 *  Created on: Sep 13, 2016
 *      Author: kotarokelley
 */
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "NumericalRecipes.h"
#include "DataStructures.h"

using namespace std;


bool VERBOSE = false;
bool PASS = false;
int passes;
int fails;

int main(int argc, char* argv[]){
	if (argc == 2){
		if (!strcmp(argv[1], "-v")|| !strcmp(argv[1], "v")|| !strcmp(argv[1], "-V")||
				!strcmp(argv[1], "V") ||!strcmp(argv[1], "verbose"))
				VERBOSE = true;
		else if (!strcmp(argv[1], "-h")|| !strcmp(argv[1], "h")|| !strcmp(argv[1], "-H")||
				!strcmp(argv[1], "H")|| !strcmp(argv[1], "help")){
				cout<<"usage: "<< argv[0] <<" < -h: help>"<< "< -v: verbose>\n>";
				exit(0);
		}
	}

	printf("****Testing NumericalRecipes Library****\n\n");
	if (VERBOSE){
		cout<<"****Testing: matrixMult()****\n";
	}
	float inA[3] = {1,2,3};
	float inB[9] = {2,1,3,3,3,2,4,1,2};
	float inC[3] = {0,0,0};
	float outC[3] = {20,10,13};
	Mat A = Mat(1,3,inA);
	Mat B = Mat(3,3,inB);
	Mat C = Mat(1,3,inC);
	if (VERBOSE){
		cout<<"A:\n";
		printMatrix(&A);
		cout<<"B:\n";
		printMatrix(&B);
		cout<<"C:\n";
		printMatrix(&C);
	}
	matrixMult(&A,&B,&C);
	if (VERBOSE){
		cout<<"A*B=C\n";
		cout<<"C after multiplication:\n";
		printMatrix(&C);
	}
	// Test if matrixMult() worked as expected.
	int counter = 0;
	for (int i=0; i<C.rows*C.cols;i++){
		if (C.dat[i]==outC[i])
			counter++;
	}

	// Test to see exception is raised when trying B*A=C.
	if (VERBOSE)
		cout<<"Testing to see if matrixMult() correctly raises exception.\n";
	int status = matrixMult(&B,&A,&C);
	if (VERBOSE){
		switch(status){
		case 0:
			cout<<"Successfully raised matrix_compatibility_exception.\n";
			break;
		case -1:
			cout<<"Failed to raise matrix_compatibility_exception. Something else went wrong.\n";
			break;
		}
	}
	if (counter==3 && status == 0){
		PASS = true;
		cout<<"\nmatrixMult(): pass.\n";
		passes++;
	}
	else{
		PASS = false;
		cout<<"\nmatrixMult(): fail.\n";
		fails++;
	}
}

