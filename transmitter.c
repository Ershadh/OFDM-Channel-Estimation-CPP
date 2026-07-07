// Replace the entire contents of transmitter.c

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include <stdbool.h>

#include "transmitter.h"
#include "config.h"

// --- Private Helper Functions ---
static void generate_qpsk_symbols(Complex* symbol_buffer, int count) {
    const double norm = 1.0 / sqrt(2.0);
    for (int i = 0; i < count; ++i) {
        int rand_val = rand() % 4;
        switch (rand_val) {
            case 0: symbol_buffer[i] = (Complex){norm, norm}; break;
            case 1: symbol_buffer[i] = (Complex){-norm, norm}; break;
            case 2: symbol_buffer[i] = (Complex){-norm, -norm}; break;
            case 3: symbol_buffer[i] = (Complex){norm, -norm}; break;
        }
    }
}

static void generate_chu_sequence(Complex* pilot_buffer, int count, int total_size) {
    for (int i = 0; i < count; ++i) {
        double phase = M_PI * i * i / (double)total_size;
        pilot_buffer[i].real = cos(phase);
        pilot_buffer[i].imag = sin(phase);
    }
}

// --- Public Functions ---

void generate_ofdm_tdm_signal(Complex* output_signal, Complex* original_symbols_per_slot) {
    fftw_complex* freq_slot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    fftw_complex* time_slot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    fftw_plan ifft_plan = fftw_plan_dft_1d(NM, freq_slot, time_slot, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int k = 0; k < K; ++k) {
        if (k == K - 1) { // Pilot slot
            generate_chu_sequence((Complex*)freq_slot, NM, NC);
        } else { // Data slots
            generate_qpsk_symbols((Complex*)freq_slot, NM);
        }
        // Store original symbols for BER calculation
        memcpy(original_symbols_per_slot + k * NM, freq_slot, sizeof(Complex) * NM);

        fftw_execute(ifft_plan);

        for (int t = 0; t < NM; ++t) {
            int frame_index = k * NM + t;
            output_signal[frame_index].real = time_slot[t][0] / (double)NM;
            output_signal[frame_index].imag = time_slot[t][1] / (double)NM;
        }
    }

    fftw_destroy_plan(ifft_plan);
    fftw_free(freq_slot);
    fftw_free(time_slot);
}

void generate_conventional_ofdm_signal(Complex* output_signal, Complex* original_symbols,
                                       int num_data_symbols, int num_pilot_symbols,
                                       const int* pilot_indices) {
    fftw_complex* freq_symbol = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* time_symbol = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_plan ifft_plan = fftw_plan_dft_1d(NC, freq_symbol, time_symbol, FFTW_BACKWARD, FFTW_ESTIMATE);

    Complex* data_symbols = (Complex*)malloc(sizeof(Complex) * num_data_symbols);
    Complex* pilot_symbols = (Complex*)malloc(sizeof(Complex) * num_pilot_symbols);
    generate_qpsk_symbols(data_symbols, num_data_symbols);
    generate_chu_sequence(pilot_symbols, num_pilot_symbols, NC);

    int data_idx = 0;
    int pilot_idx = 0;
    for (int i = 0; i < NC; ++i) {
        bool is_pilot = false;
        for (int p = 0; p < num_pilot_symbols; ++p) {
            if (i == pilot_indices[p]) {
                ((Complex*)freq_symbol)[i] = pilot_symbols[pilot_idx];
                original_symbols[i] = pilot_symbols[pilot_idx]; // Store for BER calc if needed
                pilot_idx++;
                is_pilot = true;
                break;
            }
        }
        if (!is_pilot) {
            ((Complex*)freq_symbol)[i] = data_symbols[data_idx];
            original_symbols[i] = data_symbols[data_idx];
            data_idx++;
        }
    }

    fftw_execute(ifft_plan);

    for (int t = 0; t < NC; ++t) {
        output_signal[t].real = time_symbol[t][0] / (double)NC;
        output_signal[t].imag = time_symbol[t][1] / (double)NC;
    }

    free(data_symbols);
    free(pilot_symbols);
    fftw_destroy_plan(ifft_plan);
    fftw_free(freq_symbol);
    fftw_free(time_symbol);
}

void add_cyclic_prefix(const Complex* input_signal, Complex* output_signal_with_gi, int signal_len, int gi_len) {
    memcpy(output_signal_with_gi, input_signal + signal_len - gi_len, sizeof(Complex) * gi_len);
    memcpy(output_signal_with_gi + gi_len, input_signal, sizeof(Complex) * signal_len);
}

void apply_hpa_soft_limiter(Complex* signal, int signal_length) {
    const double beta = pow(10.0, BETA_DB / 20.0);
    for (int i = 0; i < signal_length; ++i) {
        double mag = complex_mag(signal[i]);
        if (mag > beta) {
            signal[i] = complex_scale(signal[i], beta / mag);
        }
    }
}
