

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "simulation.h"
#include "config.h"
#include "complex_math.h"
#include "transmitter.h"
#include "channel.h"
#include "receiver.h"


static void qpsk_demodulate(const Complex* symbols, int* bits, int num_symbols) {
    for (int i = 0; i < num_symbols; ++i) {
        bits[2 * i] = (symbols[i].real > 0) ? 0 : 1;
        bits[2 * i + 1] = (symbols[i].imag > 0) ? 0 : 1;
    }
}

static long count_bit_errors(const int* original_bits, const int* received_bits, int num_bits) {
    long errors = 0;
    for (int i = 0; i < num_bits; ++i) {
        if (original_bits[i] != received_bits[i]) {
            errors++;
        }
    }
    return errors;
}

void run_ber_simulation(SimMode mode) {

    printf("\n--- Running Simulation For: ");
    switch(mode) {
        case MODE_IDEAL_CE: printf("Ideal Channel Estimation ---\n"); break;
        case MODE_OFDM_TDM_PROPOSED: printf("OFDM/TDM with Proposed CE (Paper's Method) ---\n"); break;
        case MODE_OFDM_FDM_PILOT: printf("Conventional OFDM with FDM-Pilot ---\n"); break;
        case MODE_OFDM_TDM_PILOT: printf("Conventional OFDM with TDM-Pilot ---\n"); break;
    }
    printf("Eb/No (dB)\tBER\n");
    printf("---------------------------\n");


    int num_pilots = NM;
    int pilot_indices[num_pilots];
    for(int i = 0; i < num_pilots; ++i) pilot_indices[i] = i * K;
    int num_data_syms = NC - num_pilots;
    int num_data_bits = num_data_syms * 2; 


    for (double eb_no_db = EB_NO_START_DB; eb_no_db <= EB_NO_END_DB; eb_no_db += EB_NO_STEP_DB) {
        long total_bits = 0;
        long total_errors = 0;

        OTDM_ReceiverState* otdm_state = NULL;
        if (mode == MODE_OFDM_TDM_PROPOSED) {
            otdm_state = otdm_receiver_init();
        }

        for (int frame = 0; frame < NUM_FRAMES_PER_TRIAL; ++frame) {

            Complex channel_ir[L];
            Complex H_est[NC];
            double noise_power_est = 0.0;

            Complex* tx_signal_base = (Complex*)malloc(sizeof(Complex) * NC);
            Complex* original_symbols = (Complex*)malloc(sizeof(Complex) * NC);
            
            if (mode == MODE_OFDM_TDM_PROPOSED) {
                generate_ofdm_tdm_signal(tx_signal_base, original_symbols);
            } else { 
                generate_conventional_ofdm_signal(tx_signal_base, original_symbols, num_data_syms, num_pilots, pilot_indices);
            }
            
            Complex* tx_signal_with_gi = (Complex*)malloc(sizeof(Complex) * (NC + GI_LENGTH));
            add_cyclic_prefix(tx_signal_base, tx_signal_with_gi, NC, GI_LENGTH);
            apply_hpa_soft_limiter(tx_signal_with_gi, NC + GI_LENGTH);

   
            Complex* rx_signal = (Complex*)malloc(sizeof(Complex) * (NC + GI_LENGTH));
            generate_rayleigh_channel_ir(channel_ir);
            apply_channel_and_noise(tx_signal_with_gi, NC + GI_LENGTH, rx_signal, channel_ir, eb_no_db);

   
            if (mode == MODE_IDEAL_CE) {
                estimate_channel_ideal(H_est, channel_ir);
                double eb_no_lin = pow(10.0, eb_no_db / 10.0);
                noise_power_est = 1.0 / eb_no_lin; 
            } else if (mode == MODE_OFDM_TDM_PROPOSED) {
                estimate_channel_ofdm_tdm(H_est, &noise_power_est, otdm_state, rx_signal);
            } else if (mode == MODE_OFDM_FDM_PILOT) {

                fftw_complex* rx_time_no_gi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
                fftw_complex* rx_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
                fftw_plan plan = fftw_plan_dft_1d(NC, rx_time_no_gi, rx_freq, FFTW_FORWARD, FFTW_ESTIMATE);
                memcpy(rx_time_no_gi, rx_signal + GI_LENGTH, sizeof(Complex) * NC);
                fftw_execute(plan);
                estimate_channel_fdm_pilot(H_est, &noise_power_est, (Complex*)rx_freq, num_pilots, pilot_indices);
                fftw_destroy_plan(plan);
                fftw_free(rx_time_no_gi);
                fftw_free(rx_freq);
            } else if (mode == MODE_OFDM_TDM_PILOT) {

                estimate_channel_ideal(H_est, channel_ir);

                double eb_no_lin = pow(10.0, eb_no_db / 10.0);
                double n0 = 1.0 / eb_no_lin;
                noise_power_est = n0;
            }
            
            Complex* equalized_symbols = (Complex*)malloc(sizeof(Complex) * NC);
            if(mode == MODE_OFDM_TDM_PROPOSED){
                perform_mmse_fde_proposed(equalized_symbols, rx_signal, H_est, noise_power_est);
            } else {            
            perform_mmse_fde(equalized_symbols, rx_signal, H_est, noise_power_est);
            }

            if (mode == MODE_OFDM_TDM_PROPOSED) {
                int num_data_syms_otdm = (K - 1) * NM;
                int num_data_bits_otdm = num_data_syms_otdm * 2;

                int* original_bits = (int*)malloc(sizeof(int) * num_data_bits_otdm);
                int* received_bits = (int*)malloc(sizeof(int) * num_data_bits_otdm);
                Complex* original_data_syms = (Complex*)malloc(sizeof(Complex) * num_data_syms_otdm);
                Complex* received_data_syms = (Complex*)malloc(sizeof(Complex) * num_data_syms_otdm);


                memcpy(original_data_syms, original_symbols, sizeof(Complex) * num_data_syms_otdm);
                memcpy(received_data_syms, equalized_symbols, sizeof(Complex) * num_data_syms_otdm);

                qpsk_demodulate(original_data_syms, original_bits, num_data_syms_otdm);
                qpsk_demodulate(received_data_syms, received_bits, num_data_syms_otdm);

                total_bits += num_data_bits_otdm;
                total_errors += count_bit_errors(original_bits, received_bits, num_data_bits_otdm);

                free(original_bits);
                free(received_bits);
                free(original_data_syms);
                free(received_data_syms);

            } else { 
                int* original_bits = (int*)malloc(sizeof(int) * num_data_bits);
                int* received_bits = (int*)malloc(sizeof(int) * num_data_bits);
                Complex* original_data_syms = (Complex*)malloc(sizeof(Complex) * num_data_syms);
                Complex* received_data_syms = (Complex*)malloc(sizeof(Complex) * num_data_syms);

                int data_idx = 0;
                for (int i = 0; i < NC; ++i) {
                    bool is_pilot = false;
                    for (int p = 0; p < num_pilots; ++p) {
                        if (i == pilot_indices[p]) {
                            is_pilot = true;
                            break;
                        }
                    }

                    if (!is_pilot) {
                        if (data_idx < num_data_syms) {
                            original_data_syms[data_idx] = original_symbols[i];
                            received_data_syms[data_idx] = equalized_symbols[i];
                            data_idx++;
                        }
                    }
                }

                qpsk_demodulate(original_data_syms, original_bits, num_data_syms);
                qpsk_demodulate(received_data_syms, received_bits, num_data_syms);

                total_bits += num_data_bits;
                total_errors += count_bit_errors(original_bits, received_bits, num_data_bits);

                free(original_bits);
                free(received_bits);
                free(original_data_syms);
                free(received_data_syms);
            }

            free(tx_signal_base); free(original_symbols); free(tx_signal_with_gi);
            free(rx_signal); free(equalized_symbols);

        }

        double ber = (double)total_errors / total_bits;
        printf("%.2f\t\t%e\n", eb_no_db, ber);

        if (otdm_state) {
            otdm_receiver_destroy(otdm_state);
        }
    }
}
