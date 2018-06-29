#include "filters.h"

//Returns sign of a number
double signum(double number)
{
	return ((number > 0) - (number < 0));
}

//Change the (double) number to an (int) written on (int) bits, knowing the (double) max allowed
double double2int(double number, double max, int bits)
{
	return signum(number)*round( (fabs(number)/max) * (pow(2,bits-1)-1) ) * max / (pow(2,bits-1)-1);
}

//Changes doubles to ints in the internal format of kissFTT
kiss_fft_cpx* doubles2ints(kiss_fft_cpx* x, int size, int bits)
{
    double max = 0.0;
    double *outR = (double*) malloc(size*sizeof(double));
    double *outRi = (double*) malloc(size*sizeof(double));

    for (int i = 0; i < size; ++i)
    {
        if (fabs(x[i].r)> max)
            max = fabs(x[i].r);
        //The following line should not change anything, since our numbers are real
        if (fabs(x[i].i)> max)
            max = fabs(x[i].i);
    }
    
    for (int i = 0; i < size; ++i)
    {
        outR[i]  = double2int(x[i].r, max, bits);
        outRi[i] = double2int(x[i].i, max, bits);
    }
    
    kiss_fft_cpx* ret;
    ret = copycpxi(outR, outRi, size);
    
    return ret;
}

kiss_fft_cpx* copycpx(double *mat, int nframe)
{
    int i;
    kiss_fft_cpx* mat2;
    mat2=(kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*nframe);
    kiss_fft_scalar zero;
    memset(&zero, 0, sizeof(zero));
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
    kiss_fft_cpx* mat2;
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

/// CALCULATING THE COEFFS FROM TRANSMITANCE VALUES IN KISSFTT STRUCT
kiss_fft_cpx* findCoeffsFIR(double* transmitance, int transmitanceLength)
{
	
    // Allocate memory to hold input and output data
    kiss_fft_cpx *coeffs, *cpx_buf;
    coeffs = (kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*transmitanceLength);
    
    // Forward transform conf
    kiss_fft_cfg cfg = kiss_fft_alloc(transmitanceLength,0,NULL,NULL);
    cpx_buf = copycpx(transmitance, transmitanceLength);
    
    // Perform fft
    kiss_fft(cfg, cpx_buf, coeffs);
    
    //Release
    free(cfg);
    free(cpx_buf);
    return coeffs;
}

/// CALCULATING THE TRANSMITANCE FROM COEFFS IN KISSFTT STRUCT
kiss_fft_cpx* transFromCoeffs(kiss_fft_cpx *coeffs, int transmitanceLength)
{
	
    // Allocate memory to hold input and output data
    kiss_fft_cpx *trans, *cpx_buf;
    trans = (kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*transmitanceLength);
    
    // Forward transform conf
    kiss_fft_cfg cfg = kiss_fft_alloc(transmitanceLength,1,NULL,NULL);

    // Perform fft
    kiss_fft(cfg, coeffs, trans);
    
    //Release
    free(cfg);
    return trans;
}

