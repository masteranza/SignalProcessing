#include "filters.h"

// Simply input number of bits and transmitance values as parameters and it'll create a symmetric filter returning real coeffs
int main(int argc, char *argv[])
{   
    int bits = atof(argv[1]);
    int transmitanceLength = (argc-2)*2-1;
    printf("Creating filter of length %d and %d bits \n", transmitanceLength, bits);
    double *transmitance = (double*) malloc(transmitanceLength*sizeof(double));
    for (int i = 0; i < argc-2; ++i)
    {
        transmitance[i] = atof(argv[i+2]);
    }
    for (int i = 0; i < argc-3; i++)
    {
        transmitance[transmitanceLength/2 + i+1] = transmitance[transmitanceLength/2 - i ];
    }
    kiss_fft_cpx *coeffs = findCoeffsFIR(transmitance, transmitanceLength);
    kiss_fft_cpx *coeffsInt = doubles2ints(coeffs, transmitanceLength, bits);
    for (int i = 0; i < transmitanceLength; ++i)
    {
        printf("%f %f\n", coeffs[i].r, coeffs[i].i);
        printf("%f %f\n", coeffsInt[i].r, coeffsInt[i].i);
    }
    return 0;
}