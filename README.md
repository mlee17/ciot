# ciot
 contrast-invariant orientation tuning

## experimental design
Stimulus patches consisting of two superimposed gratings each modulated at 3Hz (**F1**) and 5Hz (**F2**). The orientation of the **F2** stimulus changed relatively from the fixed **F1** stimulus orientation in different conditions. Two intermodulation (**IM**) terms ***F1+F2*=8Hz** and ***F2-F1*=2Hz** are expected, and the **IM** terms are expected to be greatest in the zero orientation offset condition. 

Three different contrast levels were used, and there were also baseline conditions where only F1 or F2 stimuli were presented.
10 trials per each condition, and each trial was 10 sec long.

## dataset
27 conditions in total ( (7 orientation offsets + 2 baseline conditions) * 3 contrast levels )  

Cond 1-9: High Contrast (30%)  
Cond 10-18: Med Contrast (15%)  
Cond 19-27: Low Contrast (7.5%)  

Cond 1: 0 deg offset (Highest IM response expected)  
Cond 2: 7 deg  
Cond 3: 15 deg  
Cond 4: 30 deg  
Cond 5: 45 deg  
Cond 6: 75 deg  
Cond 7: 90 deg  
Cond 8: F1(3Hz) only  
Cond 9: F2(5Hz) only  

### .mat files
**1. Axx data (frequency domain)**  
  - **Axx_c0XX_trials.mat** contains trial data for Condition **XX**.  
  -- ***Amp*** contains fourier decomposed amplitudes for a frequency range from 0Hz to 50Hz (0.5Hz step).  
  -- ***Cos*** and ***Sin*** contains the real and imaginary components from the fourier analysis, respectively.  
  -- The dimension of the data for *Amp, Cos, Sin* is 105 X 128 X 10 (Frequencies X Channels X Trials).  
  - **Axx_c0XX.mat** contains weighted averaged data across trials for Condition **XX**.  
  -- The dimension of the data for *Amp, Cos, Sin* is 105 X 128 (Frequencies X Channels).  

**2. Raw data (time domain)**  
  - **Raw_c0XX_t0YY.mat** contains data for trial number YY from Condition XX.  
