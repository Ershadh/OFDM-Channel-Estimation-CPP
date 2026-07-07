#ifndef CHANNEL_H
#define CHANNEL_H

#include "complex_math.h"

void generate_rayleigh_channel_ir(Complex* h_ir);

void apply_channel_and_noise(const Complex* tx_signal, int tx_len, Complex* rx_signal, const Complex* h_ir, double eb_no_db);

#endif 
