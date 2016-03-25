#include <stdio.h>
#include <math.h>

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
    double t = 0.123456789123456789;
    Fixed f = FLOAT2FIXED(t);
    double t2 = FIXED2DOUBLE(f);

    printf("%f %f %d", t, t2, f);
    return 0;
}