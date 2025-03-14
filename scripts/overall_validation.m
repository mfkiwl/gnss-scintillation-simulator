%% Initialization
clear all; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Model Setup
general_params = get_general_parameters();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);
irr_params_set = get_irregularity_parameters();

extrapolated_irr_params = struct('Strong', [], 'Moderate', [], 'Weak', []);

% Note: The `rhof_veff_ratio_vector` array do not depend on the
% irregularity parameters U, p1, p2 and mu0. Therefore, it is sufficient to
% extrapolate it only for one scintillation scenario, which is in this case
% the Strong case.
[extrapolated_irr_params.Strong, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params_set.Strong, general_params, rhof_veff_ratio_L1);
[extrapolated_irr_params.Moderate, ~] = ...
    freq_extrapolate(irr_params_set.Moderate, general_params, rhof_veff_ratio_L1);
[extrapolated_irr_params.Weak, ~] = ...
    freq_extrapolate(irr_params_set.Weak, general_params, rhof_veff_ratio_L1);

%% Scintillation Time Series Generation
scenarios = fieldnames(extrapolated_irr_params);
frequencies = {'L1', 'L2', 'L5'};

seed = 1;

scint_field_struct = struct();
norm_phase_sdf_struct = struct();
detrended_phase_realization_struct = struct();
mu_struct = struct();
doppler_frequency_struct = struct();
% Generate and store scintillation time series
for i = 1:numel(scenarios)
    scenario = scenarios{i};
    for j = 1:numel(frequencies)
        freq = frequencies{j};
        [scint_field_struct.(scenario).(freq), ...
            norm_phase_sdf_struct.(scenario).(freq), ...
            detrended_phase_realization_struct.(scenario).(freq), ...
            mu_struct.(scenario).(freq), ...
            doppler_frequency_struct.(scenario).(freq)] = ...
        get_scintillation_time_series( ...
            general_params, ...
            extrapolated_irr_params.(scenario).(freq), ...
            rhof_veff_ratio_vector(j), ...
            seed);
    end
end

%% Plotting functions for validation
% Time vector at the using the simulation_time parameter.
% Note: The scint_field_struct data have a different size than time_vector,
% given that its size is estimated using the helping function nicefftnum.m.
time_vector = 0: general_params.dt : general_params.simulation_time - general_params.dt;

plot_all_amp_phase_sdfs(scint_field_struct, extrapolated_irr_params, detrended_phase_realization_struct, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector);
plot_all_magnitude_phase_time_series(scint_field_struct, time_vector);