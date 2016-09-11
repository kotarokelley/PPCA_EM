/*
 * Emdata created by Kotaro Kelley 150901
 * implements 2D and 3D EM data classes
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <fftw3.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "EMdata.h"



//#include <stdint.h>
//Image::Image(){};

Image::Image(float* f_data, int* f_dim, int f_numPixels, float f_pixelWidth){
	this->rows = f_dim[0];
	this->cols = f_dim[1];
	this->depth = f_dim[2];
	this->numPixels = f_numPixels;
	
	this->data = f_data;
	
	this->update();
	
	this->origin[0] = this->rows/2;		// set origin to the the center of image
	this->origin[1] = this->cols/2;
	this->origin[2] = this->depth/2;
}

Image::Image(float* f_data, int* f_dim, int f_numPixels, int* origin, float f_pixelWidth){
	Image(f_data,f_dim, f_numPixels, f_pixelWidth);
	this->origin[0] = origin[0];
	this->origin[1] = origin[1];
	this->origin[2] = origin[2];
}

Image::Image(char* f_data, int* f_dim, int f_numPixels, float f_pixelWidth, int f_offset){
	this->rows = f_dim[0];
	this->cols = f_dim[1];
	this->depth = f_dim[2];
	this->numPixels = f_numPixels;
	this->data = new float[this->numPixels];
	std::memcpy(this->data,f_data+f_offset,(this->numPixels)*sizeof(float));
	
	this->update();
	
	this->origin[0] = this->rows/2;		// set origin to the the center of image
	this->origin[1] = this->cols/2;
	this->origin[2] = this->depth/2;	
}

Image::Image(char* f_data, int* f_dim, int f_numPixels, int* origin, float f_pixelWidth, int f_offset){
	Image(f_data,f_dim, f_numPixels, f_pixelWidth, f_offset);
	this->origin[0] = origin[0];
	this->origin[1] = origin[1];
	this->origin[2] = origin[2];
}

void Image::display(void){
	cv::Mat ogImg(cv::Size(this->cols,this->rows),CV_32FC1,this->data);
	cv::Mat reImg;
	cv::resize(ogImg,reImg,cv::Size(this->cols*4,this->rows*4));
	double min;
	double max;
	cv::minMaxLoc(reImg,&min,&max);
	cv::Mat adjMap;
	reImg.convertTo(adjMap,CV_8U,255.0/(max-min),-min*225.0/(max-min));
	//cv::convertScaleAbs(reImg,adjMap,255/max);
	
	cv::namedWindow("image",cv::WINDOW_NORMAL);
	cv::imshow("image",adjMap);
	cv::waitKey(0);
}

float* Image::getData(void){
	return data;
}

int Image::getCols(void){
	return cols;
}

int Image::getRows(void){
	return rows;
}

int Image::getDepth(void){
	return depth;
}

void Image::calcMean(void){
	float temp = 0;
	for (int i=0;i<numPixels;i++){
		temp += this->data[i];
	}
	temp /= numPixels;
	this->mean = temp;
}

float Image::getMean(void){
	return mean;
}

void Image::calcStd(void){
	float temp = 0;
	for (int i=0;i<numPixels;i++){
		temp += pow(this->data[i],2);
	}
	temp /= numPixels;
	temp -= pow(this->mean,2);
	temp = sqrt(temp);
	this->std = temp;
}

float Image::getStd(void){
	return std;
}

float Image::getPixelWidth(void){
	return pixelWidth;
}

int Image::getNumPixels(void){
	return numPixels;
}


int* Image::getOrigin(void){
	return origin;
}

void Image::changeOrigin(int newOrigin[3]){
	origin[0] = newOrigin[0];
	origin[1] = newOrigin[1];
	origin[2] = newOrigin[2];
}

void Image::absIndex( int relIndex[3], int absIndex[3]){
	//int *absIndex = new int[2];				// caller must delete
	absIndex[0] = relIndex[0] + origin[0];
	absIndex[1] = relIndex[1] + origin[1];
	//return absIndex; 
}

void Image::relIndex(int absIndex[3], int relIndex[3]){
	//int *relIndex = new int[2];				// caller must delete
	relIndex[0] = absIndex[0] - origin[0];
	relIndex[1] = absIndex[1] - origin[1];
	//return relIndex;
}

void Image::update(void){
	calcMean();
	calcStd();
	calcMin();
	calcMax();
}

void Image::calcMin(void){
	float current = data[0];
	this->minIdx = 0;
	for (int i=1; i<numPixels; i++){
		if (data[i]<current){
			current = data[i];
			this->minIdx = i;
		}
	}
	this->min = current;
}

void Image:: calcMax(void){
	float current = data[0];
	this->maxIdx = 0;
	
	for (int i=1; i<numPixels; i++){
		if (data[i]>current){
			current = data[i];
			this->maxIdx = i;
		}
	}
	this->max = current;
}

float Image:: getMin(void){
	return this->min;
}

float Image:: getMax(void){
	return this->max;
}
int	Image::getMinIdx(void){
	return this->minIdx;
}

int	Image::getMaxIdx(void){
	return this->maxIdx;
}

void Image::zeroNorm(void){
	for (int i=0; i<numPixels; i++){
		this->data[i] -= this->mean;
		data[i] /= this->std;
	}
	this->update();
}

void Image::zeroFloor(void){
	//if (this->min < 0){	
		for (int i=0; i<numPixels; i++){
			data[i] -= this->min;
		}
		this->update();
	//}
}

void Image::reScale(float newMax, float newMin){
	for (int i=0; i<this->numPixels; i++){
		this->data[i] = (this->data[i]-this->min)*((newMax-newMin)/
				(this->max-this->min))+newMin;
	}
	this->update();
}

void Image::printInfo(void){

	printf("This image has %i ",numPixels);
	printf("pixels.\n");
	printf("The first pixel has a value of %e\n",data[0]);
	printf("The second pixel has a value of %e\n",data[1]);
	printf("The min pixel value is: %e ", this->getMin());
	printf("at index: %i\n", this->getMinIdx());
	printf("The max pixel value is: %e ", this->getMax());
	printf("at index: %i\n", this->getMaxIdx());
	printf("The mean pixel value is: %e\n",this->getMean());
	printf("The std pixel value is: %e\n", this->getStd());
}
void Image::setData(float* inPointer){
	delete [] this->data;		// dealocate previous data address
	this->data = inPointer;		// assign data pointer to inPointer
}
