/*
 * parser.h created by Kotaro Kelley 150820
 * implements a parser class with several subclasses to read various input files
 */
#ifndef PARSER_H
#define PARSER_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>


class Parser{
	public:
		virtual ~ Parser(){

		}
		/*---get entire header---*/
		virtual char * getHdr(void) = 0;
		/*---parse the header---*/
		virtual void parseHdr(void) = 0;
		/*---display to screen header info---*/
		virtual void hdrInfo(void) = 0;
		/*---getData the file--should be main function to call*/
		
		virtual char * getData(bool verbose = false) = 0;
		
		virtual std::vector<int> getDim(void) = 0;
		
		virtual int32_t getPixelType(void) = 0;
		
		virtual std::vector<float> getPixelWidth(void) = 0;
		
	public:
		char  					* fileName;
		FILE  					* Input;
		unsigned long			fileSize;
		char					* data;
};

class mrcParser: public Parser{
	
public:
	~ mrcParser(){
		if (this->bmpInput!=NULL)
			fclose(this->bmpInput);
	}
	/*---parse mrc format binary files---*/	
	mrcParser(char * inputFile, int read=1);	// default read file. Otherwise, write.
	/*---return pointer to header of file, which is assumed to be 1024 bytes long---*/
	char * getHdr(void);
	
	void parseHdr(void);
	
	void hdrInfo(void);
	
	char * getData(bool verbose=false);
		
	std::vector<int> getDim(void);
	
	int getPixelType(void);	
	
	std::vector<float> getPixelWidth(void);
	
	void writeHdr();

	void writeData(float * in_data, std::vector<int> in_dim);

public:
	
	FILE  					* bmpInput;
	int32_t 				num[3];
	int32_t					pixelType;
	int32_t					mst[3];
	int32_t					m[3];
	float					d[3];
	float					angle[3];
	int32_t					axis[3];
	float					mmm1[3];
	int16_t					type;
	int16_t					nspg;
	int32_t					next;
	int16_t					dvid;
	int16_t					nblank;
	int32_t					ntst;
	char					extra[24];   // should be 1byte unsigned int
	int16_t					NumIntegers;
	int16_t					NumFloats;
	int16_t					sub;
	int16_t					zfac;
	float					mm2[2];
	float					mm3[2];
	float					mm4[2];
	int16_t					ImageType;
	int16_t					LensNum;
	int16_t					n1;
	int16_t					n2;
	int16_t					v1;
	int16_t					v2;
	float					mm5[2];
	int16_t					NumTimes;
	int16_t					ImgSequence;
	float					tilt[3];
	int16_t					NumWaves;
	int16_t					wave[5];
	float					zxy0[3];
	int32_t					NumTitles;	
};
#endif
