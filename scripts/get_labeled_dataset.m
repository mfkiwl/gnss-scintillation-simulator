%% Initialization
clearvars; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Parameter Setup
% Amount of drift velocities linearly sampled from 25 to 125m/s
drift_velocities_amount = 3;
% Amount of carrier-to-noise ratio values varying from 25 to 45 dB-Hz
carrier_to_noise_ratios_amount = 3;

% Mask angle for the line-of-sight satellites
mask_angle_deg = 15;

data_set_params = get_dataset_params(drift_velocities_amount, carrier_to_noise_ratios_amount, mask_angle_deg);

