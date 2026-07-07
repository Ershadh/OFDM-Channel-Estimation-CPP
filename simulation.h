#ifndef SIMULATION_H
#define SIMULATION_H

typedef enum {
    MODE_IDEAL_CE,
    MODE_OFDM_TDM_PROPOSED,
    MODE_OFDM_FDM_PILOT,
    MODE_OFDM_TDM_PILOT
} SimMode;


void run_ber_simulation(SimMode mode);

#endif 
