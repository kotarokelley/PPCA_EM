/*
 * Image.h created by Kotaro Kelley 150820
 * implements an EM data class with several subclasses to read various input files
 */
#ifndef EMDATA_H
#define EMDATA_H

class EMdata {
	public:
		virtual ~ EMdata(){
			delete [] data; 	// clean up dynamically allocated data array
			data = NULL;
		}
		/*---DISPLAY TO SCREEN---*/
		virtual void 	display(void) = 0;
		
		virtual float*	getData(void) = 0;

		virtual int 	getCols(void) = 0;
		
		virtual int 	getRows(void) = 0;
		
		virtual int 	getDepth(void) = 0;
		
		virtual void 	calcMean(void) = 0;
		
		virtual float 	getMean(void) = 0;
		
		virtual void 	calcStd(void) = 0;
		
		virtual float 	getStd(void) = 0;
		
		virtual float 	getPixelWidth(void) = 0;
		
		virtual int		getNumPixels(void) = 0;

		
		virtual int* 	getOrigin(void) = 0;
		
		virtual void 	changeOrigin(int *) = 0;
		
		virtual void 	absIndex(int *, int *) = 0;
		
		virtual void 	relIndex(int *, int *) = 0;
		
		virtual void	update(void) = 0;
		
		virtual void 	calcMin(void) = 0;
		
		virtual void 	calcMax(void) = 0;
		
		virtual	float 	getMin(void) = 0;
		
		virtual	float 	getMax(void) = 0;
		
		virtual int		getMinIdx(void) = 0;
		
		virtual int	 	getMaxIdx(void) = 0;
		
		virtual void	zeroNorm(void) = 0;
		
		virtual void 	zeroFloor(void) = 0;
		
		virtual void 	reScale(float, float) = 0 ;
		
		virtual void 	printInfo(void) = 0;
		
		virtual void 	setData(float*) = 0;
		
	protected:
		int cols;
		int rows;
		int depth;
		int numPixels;
		int origin[3];
		
		float mean;
		float std;
		float min;
		float max;
		int	  minIdx;
		int	  maxIdx;
		float pixelWidth;
		/**struct multData{
			int type;
			union{
				uint8_t *			data0;
				int16_t *			data1;
				float *				data2;
				//float *				fdata;
				float _Complex * 	data4;
				//int16_t	*			s16data;
				uint16_t *			data6;
				int32_t	*			data7;
			};
		}data;
		**/
		float * data;				// other data types will be rejected
};

class Image: public EMdata{
	
	public:
		/**---Call these constructor when date has already been converted to float---**/
		Image(float* f_data, int* f_dim, int f_numPixels, float f_pixelWidth=1);
		Image(float* f_data, int* f_dim, int f_numPixels, int* origin, float f_pixelWidth=1);
		/**---Call these constructor when feeding in data that's been extracted from Mrc file ---**/
		Image(char* f_data, int* f_dim, int f_numPixels, float f_pixelWidth=1, int f_offset=0);
		Image(char* f_data, int* f_dim, int f_numPixels, int* origin, float f_pixelWidth=1, int f_offset=0);

		void 	display(void);
		float*	getData(void);
		int 	getCols(void);
		int 	getRows(void);
		int 	getDepth(void);
		void 	calcMean(void);
		float 	getMean(void);
		void 	calcStd(void);
		float 	getStd(void);
		float 	getPixelWidth(void);
		int		getNumPixels(void);
		int* 	getOrigin(void);
		void 	changeOrigin(int *);
		void 	absIndex(int *, int *);
		void 	relIndex(int *, int *);
		void 	update(void);
		void 	calcMin(void);
		void 	calcMax(void);
		float 	getMin(void);
		float 	getMax(void);
		int		getMinIdx(void);
		int 	getMaxIdx(void);
		void	zeroNorm(void);
		void 	zeroFloor(void);
		void 	reScale(float, float);
		void 	printInfo(void);
		void	setData(float*);
};





#endif
