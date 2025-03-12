%% Initialization
clearvars; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Parameter Setup
drift_velocities_amount = 2;

training_data_general_params = get_training_data_general_params(drift_velocities_amount);

mask_angle_deg = 15;
% The rhof_veff_ratio_struct should be like this:
% - fields are the cities names
% - values are a struct with the fields being the satellite PRNs and their
% - values being the calculated rhof_veff_ratio,
rhof_veff_ratios = get_all_rhof_veff_ratios(training_data_general_params, mask_angle_deg);

irr_params_set = get_irregularity_parameters();