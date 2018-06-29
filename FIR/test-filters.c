#include "filters.h"
//#define DEBUG

double* generateLowPass(int length, double fractionOfSpectrum);
double* generateHighPass(int length, double fractionOfSpectrum);
double* generateMidPass(int length, double fractionOfSpectrum, double offset);
void writeTestFiles(char* filename, double* transmitance, int transmitanceLength);
int outputResults(char* filename, double* transmitance, kiss_fft_cpx* trans, kiss_fft_cpx* transInt, kiss_fft_cpx* coeffs, kiss_fft_cpx* coeffsInt, int fftSize);

int main(int argc, char *argv[])
{   
    float fraction = 0.2;
    float offset = 0.5;
    if (argc < 2)
    {
        printf("No params used\n");
    }
    if (argc > 1)
    {
        fraction = atof(argv[1]);
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

    if (argc > 2)
    {
        filtr31 = generateMidPass(31, fraction, offset);
        filtr63 = generateMidPass(63, fraction, offset);
        filtr127 = generateMidPass(127, fraction, offset);
        filtr255 = generateMidPass(255, fraction, offset);
        filtr511 = generateMidPass(511, fraction, offset);
        filtr1023 = generateMidPass(1023, fraction, offset);
        
        writeTestFiles("filtr31m-", filtr31, 31);
        writeTestFiles("filtr63m-", filtr63, 63);
        writeTestFiles("filtr127m-", filtr127, 127);
        writeTestFiles("filtr255m-", filtr255, 255);
        writeTestFiles("filtr511m-", filtr511, 511);
        writeTestFiles("filtr1023m-", filtr1023, 1023);
    }
    else if (argc > 1)
    {
        filtr31 = generateLowPass(31, fraction);        
        filtr63 = generateLowPass(63, fraction);
		filtr127 = generateLowPass(127, fraction);
        filtr255 = generateLowPass(255, fraction);
        filtr511 = generateLowPass(511, fraction);
        filtr1023 = generateLowPass(1023, fraction);

//If we ever need to test highpass filters
// double* filtr31 = generateHighPass(31, fraction);
// double* filtr63 = generateHighPass(63, fraction);
// double* filtr127 = generateHighPass(127, fraction);
// double* filtr255 = generateHighPass(255, fraction);
// double* filtr511 = generateHighPass(511, fraction);
// double* filtr1023 = generateHighPass(1023, fraction);

        writeTestFiles("filtr31-", filtr31, 31);
        writeTestFiles("filtr63-", filtr63, 63);
        writeTestFiles("filtr127-", filtr127, 127);
        writeTestFiles("filtr255-", filtr255, 255);
        writeTestFiles("filtr511-", filtr511, 511);
        writeTestFiles("filtr1023-", filtr1023, 1023);
    }
    

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////								SAVING THE RESULTS TO FILES 							////////
///////////////////////////////////////////////////////////////////////////////////////////////////

char* filenameForBits(char* basename, int bits)
{
	char bitsString[32];
    sprintf(bitsString, "%d", bits);
    char* name= malloc(strlen(basename)+7);
    strcpy(name, basename);
    strcat(name, bitsString);
    strcat(name, ".dat");
    return name;
}

void writeTestFiles(char* filename, double* transmitance, int transmitanceLength)
{

    kiss_fft_cpx *coeffs, *coeffsInt, *trans, *transInt, *cpx_buf;
    
    coeffs = findCoeffsFIR(transmitance, transmitanceLength);
    
    for(int bits=0; bits<=64; bits+=2)
    {
    	char* name = filenameForBits(filename, bits);
        
        //Getting integer approximation of the coeff on bits
        coeffsInt = doubles2ints(coeffs, transmitanceLength, bits);
        
        trans = transFromCoeffs(coeffs, transmitanceLength);
        transInt = transFromCoeffs(coeffsInt, transmitanceLength);
    
        outputResults(name, transmitance, trans, transInt, coeffs, coeffsInt, transmitanceLength);
    }
    // Perform any cleaning up
    kiss_fft_cleanup();
}

int outputResults(char* filename, double* transmitance, kiss_fft_cpx* trans, kiss_fft_cpx* transInt, kiss_fft_cpx* coeffs, kiss_fft_cpx* coeffsInt, int fftSize)
{
    FILE* fp;
    int i;
    // Open file for writing
    fp = fopen(filename, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "outputFFT: Could open output file for writing\n");
    }
    
    double mag;
	double mag2;
    //Output results
	for (i=0 ; i<fftSize ; i++)
    {
        mag = sqrt(trans[i].r * trans[i].r + trans[i].i * trans[i].i);
        mag2 = sqrt(transInt[i].r * transInt[i].r + transInt[i].i * transInt[i].i);
        
        //1  - input transmitance
        //2  - Real part of the transmitance reconstructed from calculated coefficients
        //3  - Imaginary part of the transmitance reconstructed from calculated coefficients
		//4  - Real part of the transmitance reconstructed from calculated INTERGER coefficients
        //5  - Imaginary part of the transmitance reconstructed from calculated INTERGER coefficients
        //6  - Absolute value of calculated coefficients Sqrt(Re^2 + Im^2)
        //7  - Absolute value of calculated INTEGER coefficients Sqrt(Re^2 + Im^2)
        //8  - Index number of the coefficient
        //9  - Real part of the calculated coefficients
        //10 - Imaginary part of the calculated coefficients
        //11 - Real part of the calculated INTEGER coefficients
        //12 - Imaginary part of the calculated INTEGER coefficients
        ///////															1				2		3			4			5		6	  7	  8		9			10		11			12		
        fprintf(fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g %.17g\n", transmitance[i], trans[i].r, trans[i].i, transInt[i].r, transInt[i].i, mag, mag2, i, coeffs[i].r, coeffs[i].i, coeffsInt[i].r, coeffsInt[i].i);
    }
    fclose(fp);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///////								TRANSMITANCE GENERATING PART 							////////
///////////////////////////////////////////////////////////////////////////////////////////////////
double* generateHighPass(int length, double fractionOfSpectrum)
{
    double* filter;
    filter = (double*) malloc(length*sizeof(double));
    for (int i = 0; i < length; ++i)
    {
        filter[i] = (i*1./length>=(fractionOfSpectrum/2.)) && (i*1./length<1-(fractionOfSpectrum/2.))?1:0;
        
        #ifdef DEBUG 
        printf("%d,",(int) filter[i]);
        #endif
    }
    
    return filter;
}

double* generateMidPass(int length, double fractionOfSpectrum, double offset)
{
    double* filter;
    filter = (double*) malloc(length*sizeof(double));
    
    for (int i = 0; i < length; ++i)
    {
        if (offset< fractionOfSpectrum) printf("ERROR: Offset is smaller than fractionOfSpectrum\n");
        
        if ((i*1./length >= ((offset-fractionOfSpectrum)/2.)) && (i*1./length < (offset+fractionOfSpectrum)/2.)) //midpass
        	filter[i] = 1;
        else if ((i*1./length>=(1-(offset+fractionOfSpectrum)/2.)) && (i*1./length <(1-(offset-fractionOfSpectrum)/2.))) //symmetric part of midpass
        	filter[i] = 1;
        else 
        	filter[i] = 0;
        
        #ifdef DEBUG 
        printf("%d,",(int) filter[i]);
        #endif
    }
    
    return filter;
}

double* generateLowPass(int length, double fractionOfSpectrum)
{
    double* filter;
    filter = (double*) malloc(length*sizeof(double));
    
    int i = 0;
    for (; i < (length+1)/2; i++)
    {
        filter[i] = (i*1./length>(fractionOfSpectrum/2.))?0.0:1.0;

        #ifdef DEBUG 
        printf("%d,",(int) filter[i]);
        #endif
    }
    for (; i <= length; i++)
    {
        filter[i] = filter[length - i];

        #ifdef DEBUG 
        printf("%d,",(int) filter[i]);
        #endif
    }
    return filter;
}