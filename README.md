# OFDM-Channel-Estimation-CPP
C++ implementation of OFDM channel estimation for physical layer simulation
This repository contains a C++ implementation of an OFDM signal processing chain, specifically focusing on channel estimation. The primary focus of this initial commit is establishing the structural architecture of the OFDM signal processing chain.

Current Status: The code compiles and runs, but there are known accuracy/logic errors most likely in the noise estimation module that needs to be debugged. The noise estimation block, specifically the calculation of the noise component and the subsequent noise power variance (Equations 18 and 19 of the reference paper), currently yields NaNs.
