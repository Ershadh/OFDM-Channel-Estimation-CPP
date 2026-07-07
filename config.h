#ifndef CONFIG_H
#define CONFIG_H

#define NC 256       
#define NM 16         
#define K 16          
#define GI_LENGTH NM   

#define L 8             

#define BETA_DB 3.0     

#define GAMMA_FILTER_COEFF 0.05 

#define EB_NO_START_DB 5.0
#define EB_NO_END_DB 30.0
#define EB_NO_STEP_DB 2.5
#define NUM_FRAMES_PER_TRIAL 2000 

#if (K * NM != NC)
#error "Frame parameters K, NM, and NC are inconsistent!"
#endif

#endif 
