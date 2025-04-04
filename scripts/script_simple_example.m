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