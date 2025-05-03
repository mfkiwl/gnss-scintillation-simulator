function out = cpssm(varargin)
% cpssm Compact phase-screen-based scintillation model
%
% Description:
%   CPSSM (Compact phase-screen-based scintillation model) simulates
%   ionospheric scintillation effects on GNSS signals. It performs:
%     1. Input parsing of receiver origin, velocity, time span,
%        constellations, PRNs, simulation duration, sampling rate,
%        IPP altitude, drift velocity, and plotting options.
%     2. RINEX ephemeris loading and initial parameter setup.
%     3. Satellite scenario creation using satelliteScenario.
%     4. Receiver platform definition (static or dynamic).
%     5. LOS satellite determination and filtering by constellation and ID.
%     6. Spectral parameter extrapolation (U, μ₀) across GNSS frequencies.
%     7. Geometric computation of scintillation: time series, amplitude,
%        phase, detrended phase, normalized PSD, and effective path ratio
%        for each satellite-frequency pair.
%     8. Output assembly and optional visualization of scintillation realizations.
%
% Syntax:
%   out = cpssm()
%   out = cpssm('rx_origin', [lat lon alt], 'rx_vel', value, ...
%               'datetime', value, 'svid', value, ...
%               'sim_time', value, 't_samp', value, ...
%               'ipp_altitude', value, 'drift_velocity', value, ...
%               'is_plot', true)
% Author:
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com
%
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

%% Initialize
% Add to path
[cpssm_root_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cpssm_root_dir,'libs')));
addpath(genpath(fullfile(cpssm_root_dir,'cache')));

[parsed_argins,sim_params,log] = initialize_cpssm(cpssm_root_dir, varargin);

%% Scenario
sim_params = set_scenario(log, cpssm_root_dir, sim_params, parsed_argins);

%% Scintillation realization
out = get_scintillation(sim_params);

%% Plot output
if parsed_argins.is_plot
    plot_scintillation_psd(cpssm_root_dir, out);
    plot_scintillation_time_series(out, sim_params);
end

%% Play satellite-receiver scenarios
if parsed_argins.is_play
    play(out.satelliteScenario);
end
end