/*
 * mrcParser created by Kotaro Kelley 150820
 * implements a parser class with several subclasses to read various input files
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "Parser.h"

/*---constructor----*/
mrcParser::mrcParser(char * inputFile, int read){	// default read file. Otherwise, write.
	fileName = inputFile;						// INPUT FILE NAME
	if (read){									// Read mode.
		bmpInput = fopen(fileName, "rb");		// OPEN FILE
		if(bmpInput == NULL){
			printf("Could not open file: %s", fileName);
			exit (1);
		}
		fseek(bmpInput, 0, SEEK_END);			// SET POINTER TO END OF FILE
		fileSize = ftell (bmpInput);			// GET TOTAL FILE SIZE IN BYTES
		fseek(bmpInput, 0, SEEK_SET);		    // SET POINTER TO BEGINNING OF FILE
	}
	else{										// Write mode.
		bmpInput = fopen(fileName, "w");
		if(bmpInput == NULL){
			printf("Could not open file: %s", fileName);
			exit (1);
		}
	}
}
/*---get entire header---*/
char * mrcParser::getHdr(void){						// READ HEADER AND ASSIGN TO HEADER VARIABLE
	fseek(bmpInput, 0, SEEK_SET);
	char * header = new char[1024];
	//header = (char*) malloc(1024); 		// ALLOCATE 1024 BYTES FOR THE HEADER
	fread(header, 1024, 1, bmpInput);
	return header;
}

/*---parse the header---*/
void mrcParser::parseHdr(void){
	fseek(bmpInput, 0, SEEK_SET);
	fread(num,sizeof(num),1,bmpInput);
	fread(&pixelType,sizeof(pixelType),1,bmpInput);
	fread(mst,sizeof(mst),1,bmpInput);
	fread(m,sizeof(m),1,bmpInput);
	fread(d,sizeof(d),1,bmpInput);
	fread(angle,sizeof(angle),1,bmpInput);
	fread(axis,sizeof(axis),1,bmpInput);
	fread(mmm1,sizeof(mmm1),1,bmpInput);
	fread(&type,sizeof(type),1,bmpInput);
	fread(&nspg,sizeof(nspg),1,bmpInput);
	fread(&next,sizeof(next),1,bmpInput);
	fread(&dvid,sizeof(dvid),1,bmpInput);
	fread(&nblank,sizeof(nblank),1,bmpInput);
	fread(&ntst,sizeof(ntst),1,bmpInput);
	fread(extra,sizeof(extra),1,bmpInput);
	fread(&NumIntegers,sizeof(NumIntegers),1,bmpInput);
	fread(&NumFloats,sizeof(NumFloats),1,bmpInput);
	fread(&sub,sizeof(sub),1,bmpInput);
	fread(&zfac,sizeof(zfac),1,bmpInput);
	fread(mm2,sizeof(mm2),1,bmpInput);
	fread(mm3,sizeof(mm3),1,bmpInput);
	fread(mm4,sizeof(mm4),1,bmpInput);
	fread(&ImageType,sizeof(ImageType),1,bmpInput);
	fread(&LensNum,sizeof(LensNum),1,bmpInput);
	fread(&n1,sizeof(n1),1,bmpInput);
	fread(&n2,sizeof(n2),1,bmpInput);
	fread(&v1,sizeof(v1),1,bmpInput);
	fread(&v2,sizeof(v2),1,bmpInput);
	fread(mm5,sizeof(mm5),1,bmpInput);
	fread(&NumTimes,sizeof(NumTimes),1,bmpInput);
	fread(&ImgSequence,sizeof(ImgSequence),1,bmpInput);
	fread(tilt,sizeof(tilt),1,bmpInput);
	fread(&NumWaves,sizeof(NumWaves),1,bmpInput);
	fread(wave,sizeof(wave),1,bmpInput);
	fread(zxy0,sizeof(zxy0),1,bmpInput);
	fread(&NumTitles,sizeof(NumTitles),1,bmpInput);
}
/*---display to screen header info---*/
void mrcParser::hdrInfo(void){
	printf("width  x: %i\n",num[0]);
	printf("height y: %i\n",num[1]);
	printf("# total slices z: %i\n",num[2]);
	printf("pixelType: %i\n", pixelType);
	printf("pixel width  x (um) %f\n", d[0]);
	printf("pixel width  y (um) %f\n", d[1]);
	printf("pixel height z (um) %f\n", d[2]);
	printf("intensity min/max/mean %f",mmm1[0]); printf(" %f",mmm1[1]); printf(" %f\n",mmm1[2]);
	printf("origin (um) x/y/z     %e", zxy0[1]); printf(" %e", zxy0[2]); printf(" %e\n",zxy0[0]);
	printf("pixel data type: ");
	if (pixelType == 0)
		printf("8 bit (unsigned)\n");
	else if (pixelType == 1)
		printf("16 bit (unsigned)\n");
	else if (pixelType == 2)
		printf("32 bit (signed real)\n");
	else if (pixelType == 3)
		printf("16 bit (signed complex real)\n");
	else if (pixelType == 4)
		printf("32 bit (signed complex real)\n");
	else if (pixelType == 5)
		printf("16 bit (signed) IW_EMTOM\n");
	else if (pixelType == 6)
		printf("16 bit (unsigned short)\n");
	else if (pixelType == 7)
		printf("32 bit (signed long)\n");
	else
		printf("** undefined **\n");
}
/***---get data as a 1D array of type char---***/
char * mrcParser::getData(bool verbose){
	parseHdr();
	if (verbose){
		hdrInfo();
	}
	unsigned long dataOffset = 1024+next;
	unsigned long dataSize = fileSize - dataOffset; 	// ALLOCATE MEMORY FOR DATA
	//data = (char*) malloc(dataSize);
	char * data = new char[dataSize];				// ALLOCATE MEMORY FOR DATA BY DYNAMIC DECLARATION
	fseek(bmpInput, dataOffset, SEEK_SET);			// SHIFT POSITION TO START OF DATA
	fread(data,sizeof(char),dataSize,bmpInput);		// READ IN DATA
	//fclose(bmpInput);								// DONE GETTING DATA CLOSE OUT FILE
	return data;
}

void mrcParser::writeHdr(){

	fseek(bmpInput, 0, SEEK_SET);
	fwrite(num,sizeof(num),1,bmpInput);
	fwrite(&pixelType,sizeof(pixelType),1,bmpInput);
	fwrite(mst,sizeof(mst),1,bmpInput);
	fwrite(m,sizeof(m),1,bmpInput);
	fwrite(d,sizeof(d),1,bmpInput);
	fwrite(angle,sizeof(angle),1,bmpInput);
	fwrite(axis,sizeof(axis),1,bmpInput);
	fwrite(mmm1,sizeof(mmm1),1,bmpInput);
	fwrite(&type,sizeof(type),1,bmpInput);
	fwrite(&nspg,sizeof(nspg),1,bmpInput);
	fwrite(&next,sizeof(next),1,bmpInput);
	fwrite(&dvid,sizeof(dvid),1,bmpInput);
	fwrite(&nblank,sizeof(nblank),1,bmpInput);
	fwrite(&ntst,sizeof(ntst),1,bmpInput);
	fwrite(extra,sizeof(extra),1,bmpInput);
	fwrite(&NumIntegers,sizeof(NumIntegers),1,bmpInput);
	fwrite(&NumFloats,sizeof(NumFloats),1,bmpInput);
	fwrite(&sub,sizeof(sub),1,bmpInput);
	fwrite(&zfac,sizeof(zfac),1,bmpInput);
	fwrite(mm2,sizeof(mm2),1,bmpInput);
	fwrite(mm3,sizeof(mm3),1,bmpInput);
	fwrite(mm4,sizeof(mm4),1,bmpInput);
	fwrite(&ImageType,sizeof(ImageType),1,bmpInput);
	fwrite(&LensNum,sizeof(LensNum),1,bmpInput);
	fwrite(&n1,sizeof(n1),1,bmpInput);
	fwrite(&n2,sizeof(n2),1,bmpInput);
	fwrite(&v1,sizeof(v1),1,bmpInput);
	fwrite(&v2,sizeof(v2),1,bmpInput);
	fwrite(mm5,sizeof(mm5),1,bmpInput);
	fwrite(&NumTimes,sizeof(NumTimes),1,bmpInput);
	fwrite(&ImgSequence,sizeof(ImgSequence),1,bmpInput);
	fwrite(tilt,sizeof(tilt),1,bmpInput);
	fwrite(&NumWaves,sizeof(NumWaves),1,bmpInput);
	fwrite(wave,sizeof(wave),1,bmpInput);
	fwrite(zxy0,sizeof(zxy0),1,bmpInput);
	fwrite(&NumTitles,sizeof(NumTitles),1,bmpInput);
}
/***---write data to mrc file---***/
void mrcParser::writeData(float * in_data, int * in_dim){
	unsigned long dataOffset = 1024+next;

	this->num[0] = in_dim[0];
	this->num[1] = in_dim[1];
	this->num[2] = in_dim[2];
	this->pixelType = 2;

	fseek(bmpInput, dataOffset, SEEK_SET);			// SHIFT POSITION TO START OF DATA
	fwrite(in_data,sizeof(float), num[0]*num[1]*num[2],bmpInput);		// READ IN DATA

}

/***---get dimensions of image stack---***/
int * mrcParser::getDim(void){
	return num;
}
/***---get pixel type of image stack---***/
int mrcParser::getPixelType(void){
	return pixelType;
}
/**---get pixel size in um---*/
float * mrcParser::getPixelWidth(void){
	return d;
}

