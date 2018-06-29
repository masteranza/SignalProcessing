#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"

#define PI 3.14159265
#define SQRT2 1.41421356237

//Returns sign of a number
double signum(double number);

//Change the (double) number to an (int) written on (int) bits, knowing the (double) max allowed
double double2int(double number, double max, int bits);

//Changes doubles to ints in the internal format of kissFTT
kiss_fft_cpx* doubles2ints(kiss_fft_cpx* x, int size, int bits);

///Create a KISSFTT stucture of real numbers
kiss_fft_cpx* copycpx(double *mat, int nframe);

///Create a KISSFTT stucture of complex numbers
kiss_fft_cpx* copycpxi(double* mat, double* mati, int nframe);

///CALCULATING THE COEFFS FROM TRANSMITANCE VALUES IN KISSFTT STRUCT
kiss_fft_cpx* findCoeffsFIR(double* transmitance, int transmitanceLength);

///CALCULATING THE TRANSMITANCE FROM COEFFS IN KISSFTT STRUCT
kiss_fft_cpx* transFromCoeffs(kiss_fft_cpx *coeffs, int transmitanceLength);
