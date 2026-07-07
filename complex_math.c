#include "complex_math.h"


Complex complex_add(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}


Complex complex_mult(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

Complex complex_scale(Complex a, double s) {
    Complex result;
    result.real = a.real * s;
    result.imag = a.imag * s;
    return result;
}


double complex_mag(Complex a) {
    return sqrt(a.real * a.real + a.imag * a.imag);
}


Complex complex_conj(Complex a) {
    Complex result;
    result.real = a.real;
    result.imag = -a.imag;
    return result;
}

Complex complex_div(Complex a, Complex b) {
    Complex result;
    double denominator_mag_sq = b.real * b.real + b.imag * b.imag;

    if (denominator_mag_sq == 0) {
        result.real = 0;
        result.imag = 0;
        return result;
    }

    result.real = (a.real * b.real + a.imag * b.imag) / denominator_mag_sq;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator_mag_sq;

    return result;
}
