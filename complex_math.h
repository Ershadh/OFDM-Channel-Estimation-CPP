#ifndef COMPLEX_MATH_H
#define COMPLEX_MATH_H

#include <math.h>

typedef struct {
    double real;
    double imag;
} Complex;

Complex complex_add(Complex a, Complex b);

Complex complex_mult(Complex a, Complex b);

Complex complex_scale(Complex a, double s);

double complex_mag(Complex a);

Complex complex_conj(Complex a);

Complex complex_div(Complex a, Complex b);
#endif 
