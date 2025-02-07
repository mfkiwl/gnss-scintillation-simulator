clear all; clc;

addpath(genpath(fullfile('libs')));

%% Model Initialization
general_parameters = get_general_parameters();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_parameters);
irregularity_params = get_irregularity_parameters();
[extrapolated_irregularity_parameters, extrapolated_rhof_veff_vector] = ...
    freq_extrapolate_all(irregularity_params, general_parameters, rho_oveff_l1);

%% Scintillation Time Series Generation