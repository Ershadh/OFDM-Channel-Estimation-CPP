#ifndef RECEIVER_H
#define RECEIVER_H

#include <fftw3.h>
#include "complex_math.h"
#include "config.h"

typedef struct {
    Complex* previous_filtered_gi;
    Complex* chu_pilot_nm;
    fftw_complex* filtered_gi_time;
    fftw_complex* filtered_gi_freq;
    fftw_complex* sparse_channel_est_freq;
    fftw_complex* estimated_ir_time;
    fftw_plan plan_gi_fft;
    fftw_plan plan_channel_ifft;
} OTDM_ReceiverState;

void estimate_channel_ideal(Complex* H_est, const Complex* true_ir);

OTDM_ReceiverState* otdm_receiver_init();
void otdm_receiver_destroy(OTDM_ReceiverState* state);
void estimate_channel_ofdm_tdm(Complex* H_est, double* noise_power_est, OTDM_ReceiverState* state, const Complex* received_signal);

void estimate_channel_fdm_pilot(Complex* H_est, double* noise_power_est, const Complex* rx_freq_signal, int num_pilots, const int* pilot_indices);


void perform_mmse_fde(Complex* equalized_symbols, const Complex* received_signal, const Complex* H_est, double noise_power);
void perform_mmse_fde_proposed(Complex* equalized_symbols, const Complex* received_signal, const Complex* H_est, double noise_power);

#endif 
