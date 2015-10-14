#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"

int FFTandIFFT(double *transmitance, int transmitanceL, double sampFreq);
int outputResults(char *filename, kiss_fft_cpx *results, int resultSize, double sampFreq, int fftSize);

#include <stdio.h>
typedef int Fixed;

#define FRACT_BITS 16
#define FRACT_BITS_D2 8
#define FIXED_ONE (1 << FRACT_BITS)
#define INT2FIXED(x) ((x) << FRACT_BITS)
#define FLOAT2FIXED(x) ((int)((x) * (1 << FRACT_BITS))) 
#define FIXED2INT(x) ((x) >> FRACT_BITS)
#define FIXED2DOUBLE(x) (((double)(x)) / (1 << FRACT_BITS))
#define MULT(x, y) ( ((x) >> FRACT_BITS_D2) * ((y)>> FRACT_BITS_D2) )


int main(void)
{
	int windowLength = 201;
	double sampFreq = 44100;

	// Low and high pass filters
	double transFreq = 10000;

	// Test transmitance
	double filterChar[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};  
	
	FFTandIFFT(filterChar, 1023, sampFreq);

	

	return 0;
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

int FFTandIFFT(double *transmitance, int transmitanceL, double sampFreq)
{	
	double *in;
	//for storing fixed Point -> Double
	double *outR;
	double *outRi;
	int result = 0;

	// If the window length is short, zero padding will be used
	int fftSize = (transmitanceL < 1023) ? 1023 : transmitanceL;

	// Calculate size of result data
	int resultSize = (fftSize / 2) + 1;

	// Allocate memory to hold input and output data
	in = (double *) malloc(fftSize*sizeof(double));
	outR = (double *) malloc(fftSize*sizeof(double));
	outRi = (double *) malloc(fftSize*sizeof(double));

	// Forward transform conf
	kiss_fft_cfg fwd = kiss_fft_alloc(fftSize,0,NULL,NULL);
   
    // Backward config: double and real (fixed point)
    kiss_fft_cfg inv = kiss_fft_alloc(fftSize,1,NULL,NULL);
    kiss_fft_cfg inv2 = kiss_fft_alloc(fftSize,1,NULL,NULL);
	

	kiss_fft_cpx out[fftSize], *out2, iout[fftSize], iout2[fftSize], *cpx_buf; // out[fftSize]

	int i;
	// Copy window and add zero padding (if required)
	for (i=0 ; i<transmitanceL ; i++) in[i] = transmitance[i];
	for ( ; i<fftSize ; i++) in[i] = 0;
	
	cpx_buf = copycpx(in, fftSize);
	// Perform fft
    kiss_fft(fwd, cpx_buf, out);
    
    //Allocate our fixed point
    Fixed *filterChar=(Fixed *) malloc(resultSize*sizeof(int));
    for (i=0 ; i<resultSize ; i++)
    {
    	filterChar[i]=	FLOAT2FIXED(out[i].r); //Even though it says FLOAT it should work with doubles
    	outR[i]= FIXED2DOUBLE(filterChar[i]);
    	filterChar[i]=	FLOAT2FIXED(out[i].i);
    	outRi[i]=FIXED2DOUBLE(filterChar[i]);
    }
    out2 = copycpxi(outR, outRi, resultSize);

    kiss_fft(inv, out, iout);
    kiss_fft(inv2, out2, iout2);

    outputResults("coefs.dat", out, resultSize, sampFreq, fftSize);
    outputResults("coefsReal.dat", out2, resultSize, sampFreq, fftSize);

    outputResults("TransferF.dat", iout, resultSize, sampFreq, fftSize);
    outputResults("TransferFReal.dat", iout2, resultSize, sampFreq, fftSize);


	// Perform any cleaning up
	kiss_fft_cleanup();   
    free(fwd);
    free(inv);
    free(inv2);

	return result;
}

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
