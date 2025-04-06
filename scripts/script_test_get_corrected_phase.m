% script_test_get_corrected_phase.m
%
% Description:
% This script evaluates the functionality of a newly developed function that
% corrects the phase unwrapping provided by MATLAB's unwrap function through
% Fourier interpolation. It generates scintillation time series for the L1 frequency
% using a custom scintillation model, computes both the raw unwrapped phase and the
% corrected phase via Fourier interpolation, and identifies seeds where the phase
% correction differs from the raw phase beyond a specified tolerance.
%
% Script Sections:
%
% 1. Initialization
%    - Clears the workspace and command window.
%    - Adds the necessary library and cache paths.
%
% 2. Model Setup
%    - Loads general simulation parameters and irregularity parameters.
%    - Sets a specific irregularity parameter (U = 4 for the strong scenario).
%    - Computes the rhof/veff ratio for L1 and extrapolates irregularity parameters.
%
% 3. Scintillation Time Series Generation for L1
%    - Iterates over a range of seeds to generate scintillation time series.
%    - For each seed, computes the amplitude and raw unwrapped phase using MATLAB's unwrap.
%
% 4. Phase Correction and Evaluation
%    - Applies the custom phase correction function (get_corrected_phase) which
%      uses Fourier interpolation to refine the unwrapped phase.
%    - Optionally plots the difference between the raw and corrected phase for a
%      specific seed.
%    - Collects seeds where the maximum difference between the raw and corrected
%      phase exceeds a defined tolerance.
%
% 5. Results Display
%    - Outputs the list of seeds for which the corrected phase deviates from the raw
%      phase beyond the specified tolerance.
%
% Dependencies:
% This script relies on the following custom functions from the developed library:
%   - get_general_parameters
%   - get_irregularity_parameters
%   - get_rhof_veff_ratio
%   - freq_extrapolate
%   - get_scintillation_time_series
%   - get_corrected_phase
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

%% Initialization
clear; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));
%% Model Setup
% Load general simulation parameters and irregularity parameters
general_params = get_general_parameters();
irr_params = get_irregularity_parameters();
irr_params.strong.U = 4;
% Compute the ratio needed for extrapolation (for L1)
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);

% Extrapolate irregularity parameters for the "strong" scenario
% Note: extrapolated_strong is a struct with fields L1, L2, L5; we use L1.
[extrapolated_strong, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params.strong, general_params, rhof_veff_ratio_L1);

%% Scintillation Time Series Generation for L1
seed_stop = 1000;
seeds_diff = [];  % Initialize list of seeds with phase differences
tol = pi/2;       % Tolerance threshold for phase difference

for seed = 1:seed_stop
    [scint_field, ~, detrended_phase_realization, ~, ~] = get_scintillation_time_series( ...
        general_params, ...
        extrapolated_strong.L1, ...
        rhof_veff_ratio_vector(1), ... % L1 ratio
        seed);
    
    %% Compute amplitude and unwrapped phase from the raw scintillation field
    amplitude = abs(scint_field).';
    phase_raw = unwrap(angle(scint_field)).';
    
    %% Compute corrected phase using iterative interpolation (via interpft)
    [phase_int, n_interp] = get_corrected_phase(scint_field);
    phase_int = phase_int.';
    if seed == 157
        plot(phase_raw - phase_int);
    end
    %% Check if the corrected phase differs from the raw phase beyond tolerance
    if max(abs(phase_int - phase_raw)) > tol
        seeds_diff = [seeds_diff, seed]; %#ok<AGROW>
    end
end


%% Display the list of seeds with phase differences
if isempty(seeds_diff)
    fprintf('No seeds found where the corrected phase differs from the original phase beyond the tolerance of %g.\n', tol);
else
    fprintf('Seeds with corrected phase different from original (tolerance %g):\n', tol);
    disp(seeds_diff);
end