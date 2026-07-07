#ifndef TRANSMITTER_H
#define TRANSMITTER_H

#include "complex_math.h"

void generate_ofdm_tdm_signal(Complex* output_signal, Complex* original_symbols_per_slot);

void generate_conventional_ofdm_signal(Complex* output_signal, Complex* original_symbols,
                                       int num_data_symbols, int num_pilot_symbols,
                                       const int* pilot_indices);

void add_cyclic_prefix(const Complex* input_signal, Complex* output_signal_with_gi, int signal_len, int gi_len);
void apply_hpa_soft_limiter(Complex* signal, int signal_length);

#endif 
