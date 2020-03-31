# Information metrics estimation with Gaussian Copulas framework

- Multiples-rand/Gaussian_sim_test.m: Estimating information metrics with the framework on Gaussian Data
- Multiples-rand/Gaussian_sim_test_surrogate.m: Estimating information metrics with the framework on Gaussian Data and surrogate correction
- Multiples-rand/synth_bold_test.m: Estimating information metrics with the framework on synthetic synt_bold_sig.mat
- Multiples-rand/synth_bold_test2.m: Estimating information metrics with the framework on synthetic synt_bold_sig_bp_filtered.mat
- Multiples-rand/synth_bold_test3.m: Estimating information metrics with the framework on synthetic longer_synt_bold_sig.mat
- Multiples-rand/synth_bold_test4.m: Estimating information metrics with the framework on synthetic longer_synt_bold_sig_bp_filtered.mat

### Synthetic Data:

- synt_bold_sig.mat: 90x5990 synthetic bold signal matrix
- synt_bold_sig_bp_filtered.mat: 90x5990 synthetic filtered bold signal matrix
- longer_synt_bold_sig.mat: 90x12990 synthetic bold signal matrix
- longer_synt_bold_sig_bp_filtered.mat: 90x12990 synthetic filtered bold signal matrix

### Estimation framework

- information_metrics.m: wrapper for estimating with Gaussian Copula and JIDT infodynamics
- information_metrics_calc.m: information metrics estimation with Gaussian Copula
- information_metrics_cal_JIDT.m: information metrics estimation with JIDT infodynamics