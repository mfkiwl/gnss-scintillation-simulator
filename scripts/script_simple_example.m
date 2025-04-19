% script_simple_example.m
%
% Description:
% This script demonstrates the usage of basic functions from the custom library
% for generating synthetic ionospheric scintillation time series based on phase screen
% theory. The script performs a single iteration of the model for the 'strong' scenario
% at the L1 frequency, and it prepares amplitude and phase data for further validation
% through plotting.
%
% Script Sections:
%
% 1. Initialization
%    - Clears the workspace and command window.
%    - Adds the required library and cache paths.
%
% 2. Model Setup
%    - Retrieves general simulation parameters.
%    - Obtains the rhof/veff ratio for L1.
%    - Retrieves the irregularity parameters.
%    - Performs frequency extrapolation for the strong scintillation case.
%
% 3. Scintillation Time Series Generation - Single Iteration
%    - Selects the 'strong' scenario and L1 frequency.
%    - Generates a complex scintillation time series using a custom algorithm.
%
% 4. Plotting Preparation for Validation
%    - Constructs a time vector based on the simulation time and time step.
%    - Computes amplitude and phase from the generated scintillation field.
%
% Dependencies:
% This script relies on the following custom functions from the developed library:
%   - parse_input_args
%   - get_rhof_veff_ratio
%   - get_irregularity_parameters
%   - freq_extrapolate
%   - get_scintillation_time_series
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

%% Initialization
clearvars; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Model Setup
general_params = parse_input_args();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);
irr_params_set = get_irregularity_parameters();

extrapolated_irr_params = struct('strong', [], 'moderate', [], 'weak', []);

% Note: Extrapolate only for the Strong case (the others are computed in the full code)
[extrapolated_irr_params.strong, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params_set.strong, general_params, rhof_veff_ratio_L1);
% The moderate and weak cases are not needed for a single iteration

%% Scintillation Time Series Generation - Single Iteration
% Choose one scenario and one frequency
scenario = 'strong';
freq = 'L1';

seed = 1;

[scint_field, norm_phase_sdf, detrended_phase, mu, doppler_frequency] = ...
    get_scintillation_time_series( ...
        general_params, ...
        extrapolated_irr_params.(scenario).(freq), ...
        rhof_veff_ratio_vector(1), ...
        seed);

%% Plotting functions for validation
time_vector = 0 : general_params.dt : general_params.sim_time - general_params.dt;

amplitude = abs(scint_field.');
phase = unwrap(angle(scint_field.'));