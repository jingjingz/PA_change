# Code for manuscript "A Riemann Manifold Model Framework for Longitudinal Changes in Physical Activity Patterns"

File "fmatch_subj_local.m" contains example Matlab code for estimating the deformation between baseline and follow-up PA curves 
(parallel computing for subject #id).

File "code_data_analysis.R" contins example R code for data analysis.

File "code_sim.R" contins example code for simulations, 
including comparison with method of extracting PCs from vertical changes in PA only, as well as the two-step approach based on Wrobel (2019).

File "warp_subj_fn.R" contains a wrapper function to warp baseline to target curve using Wrobel's method, and it is called in code_sim.R.
