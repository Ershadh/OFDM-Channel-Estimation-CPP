#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "receiver.h"

static void generate_chu_sequence_local(Complex* pilot_buffer, int count, int total_size) {
    for (int i = 0; i < count; ++i) {
        double phase = M_PI * i * i / (double)total_size;
        pilot_buffer[i].real = cos(phase);
        pilot_buffer[i].imag = sin(phase);
    }
}

static void high_res_interpolate(Complex* H_est_final, const Complex* H_est_sparse, int sparse_size, int final_size) {
    fftw_complex* sparse_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sparse_size);
    fftw_complex* ir_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sparse_size);
    fftw_complex* ir_padded = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * final_size);
    fftw_complex* final_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * final_size);

    memcpy(sparse_freq, H_est_sparse, sizeof(Complex) * sparse_size);

    fftw_plan ifft_plan = fftw_plan_dft_1d(sparse_size, sparse_freq, ir_time, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan fft_plan = fftw_plan_dft_1d(final_size, ir_padded, final_freq, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(ifft_plan);

    for (int i = 0; i < sparse_size; ++i) { 
        ir_time[i][0] /= (double)sparse_size;
        ir_time[i][1] /= (double)sparse_size;
    }

    memset(ir_padded, 0, sizeof(fftw_complex) * final_size);
    memcpy(ir_padded, ir_time, sizeof(fftw_complex) * sparse_size);

    fftw_execute(fft_plan);
  

    memcpy(H_est_final, final_freq, sizeof(Complex) * final_size);

    fftw_destroy_plan(ifft_plan);
    fftw_destroy_plan(fft_plan);
    fftw_free(sparse_freq);
    fftw_free(ir_time);
    fftw_free(ir_padded);
    fftw_free(final_freq);
}

// --- Channel Estimation Implementations ---

void estimate_channel_ideal(Complex* H_est, const Complex* true_ir) {
    fftw_complex* ir_padded = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* freq_resp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_plan plan = fftw_plan_dft_1d(NC, ir_padded, freq_resp, FFTW_FORWARD, FFTW_ESTIMATE);

    memset(ir_padded, 0, sizeof(fftw_complex) * NC);
    memcpy(ir_padded, true_ir, sizeof(Complex) * L);
    fftw_execute(plan);
    for (int i=0; i < NC; i++) {
        freq_resp[i][0] /= (double) NC;
        freq_resp[i][1] /= (double) NC;
    }
    memcpy(H_est, freq_resp, sizeof(Complex) * NC);

    fftw_destroy_plan(plan);
    fftw_free(ir_padded);
    fftw_free(freq_resp);
}

OTDM_ReceiverState* otdm_receiver_init() {
    OTDM_ReceiverState* state = (OTDM_ReceiverState*)malloc(sizeof(OTDM_ReceiverState));
    state->previous_filtered_gi = (Complex*)calloc(GI_LENGTH, sizeof(Complex));
    state->chu_pilot_nm = (Complex*)malloc(sizeof(Complex) * NM);
    generate_chu_sequence_local(state->chu_pilot_nm, NM, NC);
    state->filtered_gi_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    state->filtered_gi_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    state->sparse_channel_est_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    state->estimated_ir_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    state->plan_gi_fft = fftw_plan_dft_1d(NM, state->filtered_gi_time, state->filtered_gi_freq, FFTW_FORWARD, FFTW_ESTIMATE);
    state->plan_channel_ifft = fftw_plan_dft_1d(NM, state->sparse_channel_est_freq, state->estimated_ir_time, FFTW_BACKWARD, FFTW_ESTIMATE);
    return state;
}

void otdm_receiver_destroy(OTDM_ReceiverState* state) {
    if (!state) return;
    free(state->previous_filtered_gi);
    free(state->chu_pilot_nm);
    fftw_free(state->filtered_gi_time);
    fftw_free(state->filtered_gi_freq);
    fftw_free(state->sparse_channel_est_freq);
    fftw_free(state->estimated_ir_time);
    fftw_destroy_plan(state->plan_gi_fft);
    fftw_destroy_plan(state->plan_channel_ifft);
    free(state);
}

void estimate_channel_ofdm_tdm(Complex* H_est, double* noise_power_est, OTDM_ReceiverState* state, const Complex* received_signal) {
    // Time-domain filtering
    const Complex* current_gi = received_signal;
    for (int t = 0; t < GI_LENGTH; ++t) {
        Complex term1 = complex_scale(current_gi[t], GAMMA_FILTER_COEFF);
        Complex term2 = complex_scale(state->previous_filtered_gi[t], 1.0 - GAMMA_FILTER_COEFF);
        ((Complex*)state->filtered_gi_time)[t] = complex_add(term1, term2);
    }
    memcpy(state->previous_filtered_gi, state->filtered_gi_time, GI_LENGTH * sizeof(Complex));

    // FFT filtered GI
    fftw_execute(state->plan_gi_fft);
    // for (int q = 0; q < NM; ++q) {
    //     ((Complex*)state->filtered_gi_freq)[q].real /= (double)NM;
    //     ((Complex*)state->filtered_gi_freq)[q].imag /= (double)NM;
    // }

    // Reverse modulation
    for (int q = 0; q < NM; ++q) {
        ((Complex*)state->sparse_channel_est_freq)[q] = complex_div(((Complex*)state->filtered_gi_freq)[q], state->chu_pilot_nm[q]);
    }

    high_res_interpolate(H_est, (Complex*)state->sparse_channel_est_freq, NM, NC);

    double noise_power_sum = 0;
    for (int q = 0; q < NM; ++q) {
        Complex H_sparse_q = ((Complex*)state->sparse_channel_est_freq)[q];
        Complex est_rx_pilot_comp = complex_mult(H_sparse_q, state->chu_pilot_nm[q]);
        Complex R_g_q = ((Complex*)state->filtered_gi_freq)[q];
        
        Complex N_e_q = {R_g_q.real - est_rx_pilot_comp.real, R_g_q.imag - est_rx_pilot_comp.imag};
        noise_power_sum += N_e_q.real * N_e_q.real + N_e_q.imag * N_e_q.imag;
    }
    *noise_power_est = noise_power_sum / (double)NM;
}

void estimate_channel_fdm_pilot(Complex* H_est, double* noise_power_est, const Complex* rx_freq_signal, int num_pilots, const int* pilot_indices) {
    Complex* sparse_est = (Complex*)malloc(sizeof(Complex) * num_pilots);
    Complex* pilot_seq = (Complex*)malloc(sizeof(Complex) * num_pilots);
    generate_chu_sequence_local(pilot_seq, num_pilots, NC);


    for (int i = 0; i < num_pilots; ++i) {
        int idx = pilot_indices[i];
        sparse_est[i] = complex_div(rx_freq_signal[idx], pilot_seq[i]);
    }
    

    for (int i = 0; i < num_pilots - 1; ++i) {
        int start_idx = pilot_indices[i];
        int end_idx = pilot_indices[i+1];
        Complex start_val = sparse_est[i];
        Complex end_val = sparse_est[i+1];
        for (int j = start_idx; j < end_idx; ++j) {
            double factor = (double)(j - start_idx) / (end_idx - start_idx);
            H_est[j].real = start_val.real * (1.0 - factor) + end_val.real * factor;
            H_est[j].imag = start_val.imag * (1.0 - factor) + end_val.imag * factor;
        }
    }

    for(int i = 0; i < pilot_indices[0]; ++i) H_est[i] = sparse_est[0];
    for(int i = pilot_indices[num_pilots-1]; i < NC; ++i) H_est[i] = sparse_est[num_pilots-1];


    double noise_power_sum = 0;
    for(int i = 0; i < num_pilots; ++i) {
        int idx = pilot_indices[i];
        Complex est_rx_pilot = complex_mult(H_est[idx], pilot_seq[i]);
        Complex noise_val = {rx_freq_signal[idx].real - est_rx_pilot.real, rx_freq_signal[idx].imag - est_rx_pilot.imag};
        noise_power_sum += noise_val.real * noise_val.real + noise_val.imag * noise_val.imag;
    }
    *noise_power_est = noise_power_sum / num_pilots;


    free(sparse_est);
    free(pilot_seq);
}



void perform_mmse_fde(Complex* equalized_symbols, const Complex* received_signal, const Complex* H_est, double noise_power) {
    fftw_complex* rx_time_no_gi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* rx_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* eq_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_plan plan_rx_fft = fftw_plan_dft_1d(NC, rx_time_no_gi, rx_freq, FFTW_FORWARD, FFTW_ESTIMATE);

    memcpy(rx_time_no_gi, received_signal + GI_LENGTH, sizeof(Complex) * NC);
    fftw_execute(plan_rx_fft);
   
    for (int n = 0; n < NC; ++n) {
        ((Complex*)rx_freq)[n].real /= (double)NC;
        ((Complex*)rx_freq)[n].imag /= (double)NC;
    }

    double beta_lin = pow(10.0, BETA_DB / 20.0);
    double beta_sq = beta_lin * beta_lin;
    double alpha = 1.0 - exp(-beta_sq) + sqrt(M_PI * beta_sq / 2.0) * erfc(beta_lin);

    for (int n = 0; n < NC; ++n) {
        Complex H_e_n = H_est[n];
        Complex H_e_n_conj = complex_conj(H_e_n);
        double H_e_n_mag_sq = H_e_n.real * H_e_n.real + H_e_n.imag * H_e_n.imag;
        Complex denominator = {H_e_n_mag_sq * alpha * alpha + noise_power, 0};
        Complex numerator = complex_scale(H_e_n_conj, alpha);
        Complex w_n = complex_div(numerator, denominator);
        ((Complex*)eq_freq)[n] = complex_mult(w_n, ((Complex*)rx_freq)[n]);
    }

    memcpy(equalized_symbols, eq_freq, sizeof(Complex) * NC);

    fftw_destroy_plan(plan_rx_fft);
    fftw_free(rx_time_no_gi);
    fftw_free(rx_freq);
    fftw_free(eq_freq);
}

void perform_mmse_fde_proposed(Complex* equalized_symbols, const Complex* received_signal, const Complex* H_est, double noise_power) {

    fftw_complex* rx_time_no_gi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* rx_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* eq_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_complex* recovered_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);

    fftw_complex* time_slot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);
    fftw_complex* freq_slot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NM);

    fftw_plan plan_rx_fft = fftw_plan_dft_1d(NC, rx_time_no_gi, rx_freq, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan plan_eq_ifft = fftw_plan_dft_1d(NC, eq_freq, recovered_time, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_plan plan_slot_fft = fftw_plan_dft_1d(NM, time_slot, freq_slot, FFTW_FORWARD, FFTW_ESTIMATE);


    memcpy(rx_time_no_gi, received_signal + GI_LENGTH, sizeof(Complex) * NC);
    fftw_execute(plan_rx_fft);

    double beta_lin = pow(10.0, BETA_DB / 20.0);
    double beta_sq = beta_lin * beta_lin;

    double alpha = 1.0 - exp(-beta_sq) + sqrt(M_PI * beta_sq / 2.0) * erfc(beta_lin);

    for (int n = 0; n < NC; ++n) {
        Complex H_e_n = H_est[n];
        Complex H_e_n_conj = complex_conj(H_e_n);
        double H_e_n_mag_sq = H_e_n.real * H_e_n.real + H_e_n.imag * H_e_n.imag;


        Complex denominator = {H_e_n_mag_sq * alpha * alpha + noise_power, 0};
        Complex numerator = complex_scale(H_e_n_conj, alpha);
        Complex w_n = complex_div(numerator, denominator);

        ((Complex*)eq_freq)[n] = complex_mult(w_n, ((Complex*)rx_freq)[n]);
    }

    fftw_execute(plan_eq_ifft);

    for (int t = 0; t < NC; ++t) {
        ((Complex*)recovered_time)[t].real /= (double)NC;
        ((Complex*)recovered_time)[t].imag /= (double)NC;
    }


    for (int k = 0; k < K; ++k) {

        memcpy(time_slot, recovered_time + k * NM, sizeof(Complex) * NM);

        fftw_execute(plan_slot_fft);
/*        
    for (int i = 0; i < NM; ++i) {
        ((Complex*)freq_slot)[i].real /= (double)NM;
        ((Complex*)freq_slot)[i].imag /= (double)NM;
    }*/

        memcpy(equalized_symbols + k * NM, freq_slot, sizeof(Complex) * NM);
    }

    fftw_destroy_plan(plan_rx_fft);
    fftw_destroy_plan(plan_eq_ifft);
    fftw_destroy_plan(plan_slot_fft);
    fftw_free(rx_time_no_gi);
    fftw_free(rx_freq);
    fftw_free(eq_freq);
    fftw_free(recovered_time);
    fftw_free(time_slot);
    fftw_free(freq_slot);
}
