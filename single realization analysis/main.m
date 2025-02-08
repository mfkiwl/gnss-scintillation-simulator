clear all; clc;

addpath(genpath(fullfile('libs')));

%% Model Initialization
gen_params = get_general_parameters();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(gen_params);
irr_params_set = get_irregularity_parameters();

%% Scintillation Time Series Generation
seed = 1;
scint_field = get_scintillation_time_series(gen_params, irr_params_set.Severe, rhof_veff_ratio, seed);