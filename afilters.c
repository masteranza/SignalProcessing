#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"

#define PI 3.14159265
#define SQRT2 1.41421356237
int FFTandIFFT(double *transmitance, int transmitanceL, double sampFreq);
int outputResults(char *filename, kiss_fft_cpx *results, int resultSize, double sampFreq, int fftSize);
double* generateLowPass(int length, double fractionOfSpectrum);
double* generateSmoothLowPass(int length, double fractionOfSpectrum);
double* generateHighPass(int length, double fractionOfSpectrum);
double* generateMidPass(int length, double fractionOfSpectrum, double offset);
// int getCoefs(char* filename, double *transmitance, int transmitanceL, double sampFreq);
int getCoefs(char* filename, double *transmitance,double *reftransmitance, int transmitanceL, double sampFreq);
kiss_fft_cpx* copycpxi(double *mat, double *mati, int nframe);

typedef int Fixed;

#define FRACT_BITS 1
#define FRACT_BITS_D2 0
#define FIXED_ONE (1 << FRACT_BITS)
#define INT2FIXED(x) ((x) << FRACT_BITS)
// #define FLOAT2FIXED(x) ((int)((x) * (1 << FRACT_BITS))) 
#define FIXED2INT(x) ((x) >> FRACT_BITS)
// #define FIXED2DOUBLE(x) (((double)(x)) / (1 << FRACT_BITS))
#define MULT(x, y) ( ((x) >> FRACT_BITS_D2) * ((y)>> FRACT_BITS_D2) )

double fixed2double(int x, int bits) 
{
	return (((double)(x)) / (1 << bits));
}
int double2fixed(double x, int bits) 
{
	return ((int)((x) * (1 << bits)));
}




int main(int argc, char *argv[])
{
	int windowLength = 201;
	double sampFreq = 44100;

	// Test transmitance
	double filterChar[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., \
1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, \
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
//LOW PASS
// {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
// 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
// 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  

	float part = 0.2;
	float offset = 0.5;
  if (argc < 2)
    {
      printf("No params used\n");
      
    }
if (argc > 1)
{
  
  part = atof(argv[1]); 
 if (argc > 2)
{
  
  offset = atof(argv[2]); 
}
}
	double* filtr31;
	double* filtr63;
	double* filtr127;
	double* filtr255;
	double* filtr511;
	double* filtr1023;
	
	double* reffiltr31;
	double* reffiltr63;
	double* reffiltr127;
	double* reffiltr255;
	double* reffiltr511;
	double* reffiltr1023;

	if (argc > 2)
	{
		filtr31 = generateMidPass(31, part, offset);
		filtr63 = generateMidPass(63, part, offset);
		filtr127 = generateMidPass(127, part, offset);
		filtr255 = generateMidPass(255, part, offset);
		filtr511 = generateMidPass(511, part, offset);
		filtr1023 = generateMidPass(1023, part, offset);

		getCoefs("filtr31m-", filtr31, reffiltr31, 31, sampFreq);
		getCoefs("filtr63m-", filtr63, reffiltr63, 63, sampFreq);
		getCoefs("filtr127m-", filtr127, reffiltr127, 127, sampFreq);
		getCoefs("filtr255m-", filtr255, reffiltr255, 255, sampFreq);
		getCoefs("filtr511m-", filtr511, reffiltr511, 511, sampFreq);
		getCoefs("filtr1023m-", filtr1023, reffiltr1023, 1023, sampFreq);
	}	
	else if (argc > 1)
	{
		filtr31 = generateSmoothLowPass(31, part);
		reffiltr31 = generateLowPass(31, part);
		getCoefs("filtr31-", filtr31, reffiltr31, 31, sampFreq);

		filtr63 = generateSmoothLowPass(63, part);
		reffiltr63 = generateLowPass(63, part);
		getCoefs("filtr63-", filtr63, reffiltr63, 63, sampFreq);

		filtr127 = generateSmoothLowPass(127, part);
		reffiltr127 = generateLowPass(127, part);
		getCoefs("filtr127-", filtr127, reffiltr127, 127, sampFreq);

		filtr255 = generateSmoothLowPass(255, part);
		reffiltr255 = generateLowPass(255, part);
		getCoefs("filtr255-", filtr255, reffiltr255, 255, sampFreq);

		filtr511 = generateSmoothLowPass(511, part);
		reffiltr511 = generateLowPass(511, part);
		getCoefs("filtr511-", filtr511, reffiltr511, 511, sampFreq);

		filtr1023 = generateSmoothLowPass(1023, part);
		reffiltr1023 = generateLowPass(1023, part);
		getCoefs("filtr1023-", filtr1023, reffiltr1023, 1023, sampFreq);

	}

	


	
	

	// double* filtr31 = generateHighPass(31, part);
	// double* filtr63 = generateHighPass(63, part);
	// double* filtr127 = generateHighPass(127, part);
	// double* filtr255 = generateHighPass(255, part);
	// double* filtr511 = generateHighPass(511, part);
	// double* filtr1023 = generateHighPass(1023, part);
	
	// getCoefs("filtr31h-", filtr31, 31, sampFreq);
	// getCoefs("filtr63h-", filtr63, 63, sampFreq);
	// getCoefs("filtr127h-", filtr127, 127, sampFreq);
	// getCoefs("filtr255h-", filtr255, 255, sampFreq);
	// getCoefs("filtr511h-", filtr511, 511, sampFreq);
	// getCoefs("filtr1023h-", filtr1023, 1023, sampFreq);

	// FFTandIFFT(filterChar, 1023, sampFreq);
	return 0;
}

double* generateHighPass(int length, double fractionOfSpectrum)
{
	double *filter;
	filter = (double *) malloc(length*sizeof(double));
	printf("Filter of length %d: [", length);
	for (int i = 0; i < length; ++i)
	{

		filter[i] = (i*1./length>=(fractionOfSpectrum/2.)) && (i*1./length<1-(fractionOfSpectrum/2.))?1:0;
		printf("%d,",(int) filter[i]);
	}
	printf("]\n");
	return filter;
}

double* generateMidPass(int length, double fractionOfSpectrum, double offset)
{
	double *filter;
	filter = (double *) malloc(length*sizeof(double));
	printf("Filter of length %d: [", length);
	for (int i = 0; i < length; ++i)
	{
//Uzupełnij
		if (offset< fractionOfSpectrum) printf("ERROR: Offset is smaller than fractionOfSpectrum\n");

		filter[i] = ((i*1./length>=((offset-fractionOfSpectrum)/2.)) && (i*1./length < (offset+fractionOfSpectrum)/2 ))
					||
					((i*1./length>=(1-(offset+fractionOfSpectrum)/2.)) && (i*1./length <(1-(offset-fractionOfSpectrum)/2)))?1:0;

		// filter[i] = (i*1./length>=(fractionOfSpectrum/2.)) && (i*1./length<1-(fractionOfSpectrum/2.))?1:0;
		printf("%d,",(int) filter[i]);
	}
	printf("]\n");
	return filter;
}

double* generateLowPass(int length, double fractionOfSpectrum)
{
	double *filter;
	filter = (double *) malloc(length*sizeof(double));
	// printf("Filter of length %d: [", length);
	for (int i = 0; i < length; i++)
	{

		filter[i] = (i*1./length>=(fractionOfSpectrum/2.)) && (i*1./length<1-(fractionOfSpectrum/2.))?0:1;
		// printf("%d,",(int) filter[i]);
	}

	// printf("]\n");
	return filter;
}

double* generateSmoothLowPass(int length, double fractionOfSpectrum)
{
	double *filter;
	filter = (double *) malloc(length*sizeof(double));
	// printf("Filter of length %d: [", length);
	int i = 0;
	for (; i < (length+1)/2; i++)
	{
		filter[i] = (i*1./length>(fractionOfSpectrum/2.))?0.0:1.0;
		// filter[i] = (i*1./length>(fractionOfSpectrum))?0.0:1.0;
		// filter[i] = filter[i]* cos((i*1./length)*PI/(2. * fractionOfSpectrum))*cos((i*1./length)*PI/(2. * fractionOfSpectrum));

		// if (length==255) printf("for %d of %d %f\n,",i, length, filter[i]);
	}
	for (; i <= length; i++)
	{
		filter[i] = filter[length - i];
		// if (length==255) printf("for %d of %d %f\n,",i, length, filter[i]);
	}
	// printf("]\n");
	return filter;
}

kiss_fft_cpx* copycpx(double *mat, int nframe)
{
	int i;
	kiss_fft_cpx *mat2;
	mat2=(kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*nframe);
        kiss_fft_scalar zero;
        memset(&zero,0,sizeof(zero) );
	for(i=0; i<nframe ; i++)
	{
		mat2[i].r = mat[i];
		mat2[i].i = zero;
	}
	return mat2;
}

kiss_fft_cpx* copycpxi(double *mat, double *mati, int nframe)
{
	int i;
	kiss_fft_cpx *mat2;
	mat2=(kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*nframe);
    
    kiss_fft_scalar zero;
    memset(&zero,0,sizeof(zero) );

	for(i=0; i<nframe ; i++)
	{
		mat2[i].r = mat[i];
		mat2[i].i = mati[i];
	}
	return mat2;
}

kiss_fft_cpx* doubles2fixed(kiss_fft_cpx* x, int size, int bits) 
{
	double max = 0.0;
	double *outR;
	double *outRi;
	outR = (double *) malloc(size*sizeof(double));
	outRi = (double *) malloc(size*sizeof(double));

	for (int i = 0; i < size; ++i)
	{
		if (fabs(x[i].r)> max)
			max = fabs(x[i].r);
		if (fabs(x[i].i)> max)
			max = fabs(x[i].i);
	}
	
	for (int i = 0; i < size; ++i)
	{
		outR[i]  = ((x[i].r > 0) - (x[i].r < 0))*round((fabs(x[i].r)/max) * (pow(2,bits-1)-1)) * max / (pow(2,bits-1)-1);
		outRi[i] = ((x[i].i > 0) - (x[i].i < 0))*round((fabs(x[i].i)/max) * (pow(2,bits-1)-1)) * max / (pow(2,bits-1)-1);

		// if (size == 31) printf("%f %f\n", outR[i], outR[i]);
	}

	kiss_fft_cpx* ret;
	ret = copycpxi(outR, outRi, size);

	// printf("For bits %d and size %d: (%f %f), (%f %f) \n", bits, size, x[1].r, outR[1], x[2].r, outR[2]);
	return ret;
}


int getCoefs(char* filename, double *transmitance,double *reftransmitance, int transmitanceL, double sampFreq)
{
	FILE *fp;
	double *in;
	int result;
	int i;

	

	// // Output result
	// for (i=0 ; i<resultSize ; i++)
	// {	
	// 	double freq = sampFreq * i / fftSize;
	// 	double mag = sqrt(out_cpx[i].r * out_cpx[i].r + out_cpx[i].i * out_cpx[i].i);
	// 	double magdB = 20 * log10(mag);
	// 	double phase = atan2(out_cpx[i].i, out_cpx[i].r);
	// 	fprintf(fp, "%02d %f %f %f %f\n", i, freq, mag, magdB, phase);
	// }

	// // Perform any cleaning up
	// finalise:
	// kiss_fft_cleanup();   
 	//    free(fft);
	//for storing fixed Point -> Double

	double *outR;
	double *outRi;
	result = 0;

	int fftSize =  transmitanceL;

	// Calculate size of result data
	int resultSize = (fftSize / 2) + 1;

	// Allocate memory to hold input and output data
	in = (double *) malloc(fftSize*sizeof(double));
	// outR = (double *) malloc(fftSize*sizeof(double));
	// outRi = (double *) malloc(fftSize*sizeof(double));

	// Forward transform conf
	kiss_fft_cfg fwd = kiss_fft_alloc(fftSize,0,NULL,NULL);
   
    // Backward config: double and real (fixed point)
    kiss_fft_cfg inv = kiss_fft_alloc(fftSize,1,NULL,NULL);
    kiss_fft_cfg inv2 = kiss_fft_alloc(fftSize,1,NULL,NULL);
	

	kiss_fft_cpx out[fftSize], *out2, iout[fftSize], iout2[fftSize], *cpx_buf; // out[fftSize]

	// Copy window and add zero padding (if required)
	for (i=0 ; i<transmitanceL ; i++) in[i] = transmitance[i];
	for (; i<fftSize; i++) in[i] = 0;
	
	cpx_buf = copycpx(in, fftSize);
	// Perform fft
    kiss_fft(fwd, cpx_buf, out);
    
    //Allocate our fixed point
    Fixed *filterChar=(Fixed *) malloc(fftSize*sizeof(int));
    Fixed *filterChari=(Fixed *) malloc(fftSize*sizeof(int));

    for(int bits=0; bits<=64; bits+=2)
    {
 		char bits_string[32];
		sprintf(bits_string, "%d", bits);
		
		char* name= malloc(strlen(filename)+1+6);
		strcpy(name, filename); /* copy name into the new var */
		strcat(name, bits_string);
		strcat(name, ".dat");

		
    	fp = fopen(name, "w");
		if (fp == NULL) 
		{
			result = 1;
			fprintf(stderr, "outputFFT: Could open output file for writing\n");
			// goto finalise;
		}
		
	    // for (i=0 ; i<fftSize ; i++)
	    // {
	    // 	filterChar[i]=	double2fixed(out[i].r, bits); 
	    // 	outR[i]= fixed2double(filterChar[i], bits);
	    // 	filterChari[i]=	double2fixed(out[i].i, bits);
	    // 	outRi[i]= fixed2double(filterChari[i], bits);
	    // }
	    //remember fftSize > resultSize
	    // out2 = copycpxi(outR, outRi, fftSize);
	    out2 = doubles2fixed(out, fftSize, bits);

	    kiss_fft(inv, out, iout);
	    kiss_fft(inv2, out2, iout2);

	    double mag;
	    double mag2;
	    for (i=0 ; i<fftSize ; i++)
		{
			mag = sqrt(iout[i].r * iout[i].r + iout[i].i * iout[i].i);

			// printf("-- %f %f\n", iout2[i].r,iout2[i].i);
			mag2 = sqrt(iout2[i].r * iout2[i].r + iout2[i].i * iout2[i].i);
																	//1				2			3			4			5		6	  7	   8	9			10		11			12			13
			fprintf(fp, "%e %e %e %e %e %e %e %d %e %e %e %e %e\n", transmitance[i], iout[i].r, iout[i].i, iout2[i].r, iout2[i].i, mag, mag2 ,i, out[i].r, out[i].i, out2[i].r, out2[i].i, reftransmitance[i]);
		}
	}
	// Perform any cleaning up
	kiss_fft_cleanup();   
    free(fwd);
    free(inv);
    free(inv2);
    return result;
}

// int FFTandIFFT(double *transmitance, int transmitanceL, double sampFreq)
// {	
// 	double *in;
// 	//for storing fixed Point -> Double
// 	double *outR;
// 	double *outRi;
// 	int result = 0;

// 	// If the window length is short, zero padding will be used
// 	int fftSize = (transmitanceL < 1023) ? 1023 : transmitanceL;

// 	// Calculate size of result data
// 	int resultSize = (fftSize / 2) + 1;

// 	// Allocate memory to hold input and output data
// 	in = (double *) malloc(fftSize*sizeof(double));
// 	outR = (double *) malloc(fftSize*sizeof(double));
// 	outRi = (double *) malloc(fftSize*sizeof(double));

// 	// Forward transform conf
// 	kiss_fft_cfg fwd = kiss_fft_alloc(fftSize,0,NULL,NULL);
   
//     // Backward config: double and real (fixed point)
//     kiss_fft_cfg inv = kiss_fft_alloc(fftSize,1,NULL,NULL);
//     kiss_fft_cfg inv2 = kiss_fft_alloc(fftSize,1,NULL,NULL);
	

// 	kiss_fft_cpx out[fftSize], *out2, iout[fftSize], iout2[fftSize], *cpx_buf; // out[fftSize]

// 	int i;
// 	// Copy window and add zero padding (if required)
// 	for (i=0 ; i<transmitanceL ; i++) in[i] = transmitance[i];
// 	for ( ; i<fftSize ; i++) in[i] = 0;
	
// 	cpx_buf = copycpx(in, fftSize);
// 	// Perform fft
//     kiss_fft(fwd, cpx_buf, out);
    
//     //Allocate our fixed point
//     Fixed *filterChar=(Fixed *) malloc(fftSize*sizeof(int));
//     Fixed *filterChari=(Fixed *) malloc(fftSize*sizeof(int));
//     for (i=0 ; i<fftSize ; i++)
//     {
//     	filterChar[i]=	FLOAT2FIXED(out[i].r); //Even though it says FLOAT it should work with doubles
//     	outR[i]= FIXED2DOUBLE(filterChar[i]);
//     	filterChari[i]=	FLOAT2FIXED(out[i].i);
//     	outRi[i]= FIXED2DOUBLE(filterChari[i]);
//     }
//     //remember fftSize > resultSize
//     out2 = copycpxi(outR, outRi, fftSize);
//     for (i=0 ; i<fftSize ; i++)
// 	{
// 		printf("(%f %f) (%f %f) (%d,%d)\n", out[i].r,out[i].i,out2[i].r,out2[i].i, filterChar[i],filterChari[i]);
// 	}

//     kiss_fft(inv, out, iout);
//     kiss_fft(inv2, out2, iout2);

//     outputResults("coefs.dat", out, resultSize, sampFreq, fftSize);
//     outputResults("coefsReal.dat", out2, resultSize, sampFreq, fftSize);

//     outputResults("TransferF.dat", iout, resultSize, sampFreq, fftSize);
//     outputResults("TransferFReal.dat", iout2, resultSize, sampFreq, fftSize);


// 	// Perform any cleaning up
// 	kiss_fft_cleanup();   
//     free(fwd);
//     free(inv);
//     free(inv2);

// 	return result;
// }

int outputResults(char *filename, kiss_fft_cpx *results, int resultSize, double sampFreq, int fftSize)
{
	FILE *fp;
	int i;
	// Open file for writing
	fp = fopen(filename, "w");
	if (fp == NULL) 
	{
		fprintf(stderr, "outputFFT: Could open output file for writing\n");
	}

	// Output result
	for (i=0 ; i<resultSize ; i++)
	{
		double freq = sampFreq * i / fftSize;
		double mag = sqrt(results[i].r * results[i].r + results[i].i * results[i].i);
		double magdB = 20 * log10(mag);
		double phase = atan2(results[i].i, results[i].r);
		fprintf(fp, "%02d %f %f %f %f %f %f\n", i, freq, results[i].r, results[i].i, mag, magdB, phase);
	}
	return 0;
}
