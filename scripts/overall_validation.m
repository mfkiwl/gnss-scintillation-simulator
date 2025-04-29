% overall_validation.m
%
% Description:
%
% This script generates plots comparing the behavior of synthetic ionospheric scintillation 
% time series produced by a newly developed algorithm based on phase screen theory. The generated 
% plots display the time series alongside key parameters characterizing them, and they compare 
% the estimated power spectral densities with the theoretical expectations. This script serves 
% as a reference for GNSS researchers analyzing scintillation models.
% 
% Script Sections:
% 
% 1. Initialization
%    - Clears the workspace and command window.
%    - Adds the required library paths for custom functions and cached data.
% 
% 2. Model Setup
%    - Retrieves general simulation parameters.
%    - Obtains the ratio between rhof and veff for the L1 frequency.
%    - Retrieves the irregularity parameters for different scintillation scenarios.
%    - Performs frequency extrapolation of the irregularity parameters using the Strong 
%      scintillation case as reference (since the rhof_veff_ratio_vector does not depend on 
%      specific irregularity parameters).
% 
% 3. Scintillation Time Series Generation
%    - Defines the scintillation scenarios (strong, moderate, weak) and frequency bands (L1, L2, L5).
%    - Initializes structures to store the time series and associated parameters.
%    - Generates synthetic scintillation time series by iterating over each scenario and frequency,
%      applying the custom algorithm.
% 
% 4. Plotting for Validation
%    - Creates a time vector based on the simulation time and time step.
%    - Generates plots to validate the simulation by comparing the generated time series with 
%      theoretical power spectral densities and visualizing amplitude and phase behaviors.
% 
% Dependencies:
%
% This script relies on the following custom functions from the developed library:
%   - parse_input_args
%   - get_rhof_veff_ratio
%   - get_irregularity_parameters
%   - freq_extrapolate
%   - get_scintillation_time_series
%   - plot_all_amp_phase_sdfs
%   - plot_all_magnitude_phase_time_series
%   - nicefftnum
% 
% Note:
% 
% - Detailed information regarding input parameters, expected outputs, and specific assumptions 
%   is provided within each of the individual functions.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

%% Initialization
clearvars; clc;

scint_realization = cspsm('constellation', 'gps', 'log_lvl', 'INFO','frequency', ["L1", "L2"]);

%% Plotting functions for validation
% Time vector at the using the sim_time parameter.
% Note: The scint_field_struct data have a different size than time_vector,
% given that its size is estimated using the helping function nicefftnum.m.
time_vector = 0 : general_params.dt : general_params.sim_time - general_params.dt;

plot_all_amp_phase_sdfs(scint_field_struct, extrapolated_irr_params, detrended_phase_realization_struct, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector);
plot_all_magnitude_phase_time_series(scint_field_struct, time_vector);