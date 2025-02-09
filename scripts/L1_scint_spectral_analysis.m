clear all; clc;

addpath(genpath(fullfile(pwd,'..','libs')));
addpath(fullfile(pwd,'..','cache'));

% Model Setup
gen_params = get_general_parameters();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(gen_params);
irr_params_set = get_irregularity_parameters();

% Scintillation Time Series Generation for L1 GPS frequency band
seed = 42;
[scint_field,norm_phase_sdf,detrended_phase_realization,mu,doppler_frequency] = get_scintillation_time_series(gen_params, irr_params_set.Mild, rhof_veff_ratio_L1, seed);

% Plots
plot_amp_phase_sdf(irr_params_set.Mild, scint_field, detrended_phase_realization, norm_phase_sdf, doppler_frequency, mu, rhof_veff_ratio_L1);