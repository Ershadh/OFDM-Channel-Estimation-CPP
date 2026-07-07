#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "channel.h"
#include "config.h"

static double generate_gaussian_rv() {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

void generate_rayleigh_channel_ir(Complex* h_ir) {

    double norm = 1.0 / sqrt((double)L);

    for (int i = 0; i < L; ++i) {
        h_ir[i].real = norm * generate_gaussian_rv() * (1.0 / sqrt(2.0));
        h_ir[i].imag = norm * generate_gaussian_rv() * (1.0 / sqrt(2.0));
    }
}

void apply_channel_and_noise(const Complex* tx_signal, int tx_len, Complex* rx_signal, const Complex* h_ir, double eb_no_db) {
    for (int n = 0; n < tx_len; ++n) {
        rx_signal[n].real = 0;
        rx_signal[n].imag = 0;
        for (int l = 0; l < L; ++l) {
            if (n >= l) {
                // rx[n] += h[l] * tx[n-l]
                rx_signal[n] = complex_add(rx_signal[n], complex_mult(h_ir[l], tx_signal[n - l]));
            }
        }
    }


    double signal_power = 0;
    for (int i = 0; i < tx_len; i++) {
        signal_power += tx_signal[i].real * tx_signal[i].real + tx_signal[i].imag * tx_signal[i].imag;
    }
    signal_power /= tx_len;

    double eb_no_lin = pow(10.0, eb_no_db / 10.0);
    double bits_per_symbol = 2.0; 
    double symbol_energy = signal_power;
    double bit_energy = symbol_energy / bits_per_symbol;
    double n0 = bit_energy / eb_no_lin;
    double noise_std_dev = sqrt(n0 / 2.0); 

    for (int i = 0; i < tx_len; i++) {
        rx_signal[i].real += generate_gaussian_rv() * noise_std_dev;
        rx_signal[i].imag += generate_gaussian_rv() * noise_std_dev;
    }
}
