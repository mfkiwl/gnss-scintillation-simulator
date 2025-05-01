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
%% Add to path
[cspsm_root_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cspsm_root_dir,'libs')));
addpath(genpath(fullfile(cspsm_root_dir,'cache')));

%% instantiate simulation parameters
sim_params = get_cte_sim_params();

%% handle input args
[parsed_argins, log] = parse_input_args(cspsm_root_dir, ...
    sim_params.cte.all_constellations, varargin{:});

% set basic simulation parameters directly from the input arguments: IPP
% altitude, drift velocity and receiver parameters, etc.
sim_params = set_sim_params_from_argin(sim_params, parsed_argins);

%% get RINEX file
[sim_params.time_range, rinex] = get_rinex(cspsm_root_dir, parsed_argins);

%% Get line-of-sight (LOS) satellites

% create a satellite scenario
sat_scen = satelliteScenario(sim_params.time_range.start, ...
    sim_params.time_range.end, sim_params.geo_tsamp);

% instantiate satellites for each constellation from the RINEX ephemerides
all_sats = satellite(sat_scen, rinex);

% FIXME: It is not clear to me how the velocity is defined or, at least,
% maintained constant. `geoTrajectory()` seems to require a starting and
% ending position, as well as the travel duration to compute the
% trajectory. It has some ways to add velocity (or speed) information, but
% it doesn't seem straightforward to set a constant velocity. It is because
% the velocity seems to have a strong dependence on the ending position and
% simulation time: A distant path seems to make the body accelerate during
% the trajectory. In other words, it seems to have a dynamic velocity. For
% the sake of first try, the receiver is being here to be static
% SEE: https://www.mathworks.com/help/satcom/ref/geotrajectory-system-object.html
rx_traj = geoTrajectory( ...
    [parsed_argins.rx_origin; parsed_argins.rx_origin], ...
    [0 sim_params.sim_time], ...
    'SampleRate', parsed_argins.t_samp);

rx = platform(sat_scen, rx_traj, 'Receiver');

% get parameters of the satellites in LOS with the receiver
ac = access(all_sats, rx);
los_sats_params = ac.accessIntervals;

% set the simulation parameters' frequencies and constellations based on
% the LOS satellites' parameters and user input arguments
sim_params = set_constellation_freq_svid(log, sim_params, ...
    parsed_argins, los_sats_params);

% get user-filtered LOS sats IDS (filtered by constellation and SV IDS)
filtered_los_sats_params = get_filtered_los_sat_params(log, sim_params, ...
    los_sats_params);

filtered_los_sats = all_sats(ismember(all_sats.Name, ...
    filtered_los_sats_params.Source));
sim_params.sats = filtered_los_sats;

%% Compute the nongeometric-dependent outputs
% They are
% - U and μ₀: dependent only on the frequency (computed via extrapolation)
% - p₁ and p₂: independent
% - Doppler frequency: dependent only on t_samp and FFT samples
% - UTC time: independent
% Frequency support where the PSD are plotted
nfft = nicefftnum(sim_params.sim_time / sim_params.t_samp);
out.doppler_frequency_support = (-nfft/2 : nfft/2-1) / (nfft * sim_params.t_samp);
% UTC timestamps of the propagation geometry
[~, ~, out.time_utc] = states(rx, 'CoordinateFrame','geographic');
for constellation = sim_params.constellations
    % initialize a 0x0 struct array which will store geometric-related
    % outputs. Each, struct scalar in this array is related to a sat-rx
    % geometry
    out.(constellation).scenario = struct([]);
    for i = 1:numel(sim_params.freqs.(constellation).value)
        freq_name = sim_params.freqs.(constellation).name(i);
        freq_value = sim_params.freqs.(constellation).value(i);
        
        % extrapolate U, and μ₀ from `freq_ref` to `freq`
        out.(constellation).spectral.(freq_name) = ...
            extrapolate_spectral_params(sim_params, freq_value);
    end
end

%% Compute the geometric related output
% for all filtered LOS sat (of a given contellation)
for i = 1:numel(sim_params.sats)
    % get LOS sat and constellation name
    sat = sim_params.sats(i);
    sat_constellation = sat.OrbitPropagator;
    % save sat and receiver to output
    out.(sat_constellation).scenario(end+1).sat = sat;
    out.(sat_constellation).scenario(end).rx = rx;
    % for this geometry propagation, get the ρF/veff for the reference
    % frequency
    rhof_veff_ratio_ref = get_scaling_param(rx, sat, sim_params);

    % for all valid frequencies for this satellite contellation
    for j = 1:numel(sim_params.freqs.(sat_constellation).name)
        freq_name = sim_params.freqs.(sat_constellation).name(j);
        freq_value = sim_params.freqs.(sat_constellation).value(j);

        rhof_veff_ratio = extrapolate_scaling_param(sim_params, freq_value, rhof_veff_ratio_ref);
        
        % get scintillation realization
        [scint_field, norm_phase_psd, detrended_phase, mu] = ...
            get_scintillation_time_series( ...
            out.doppler_frequency_support, ...
            out.(constellation).spectral.(freq_name), ...
            rhof_veff_ratio, ...
            nfft, ...
            sim_params.seed);

        % compute amplitude and phase time series of the received
        % scintillation signal
        amplitude = abs(scint_field.');
        phase = unwrap(angle(scint_field.'));
        
        % outputs
        out.(sat_constellation).scenario(end).(freq_name).scint_field = scint_field;
        out.(sat_constellation).scenario(end).(freq_name).amplitude = amplitude;
        out.(sat_constellation).scenario(end).(freq_name).phase = phase;
        out.(sat_constellation).scenario(end).(freq_name).detrended_phase = detrended_phase;
        out.(sat_constellation).scenario(end).(freq_name).mu = mu;
        out.(sat_constellation).scenario(end).(freq_name).norm_phase_psd = norm_phase_psd;
        out.(sat_constellation).scenario(end).(freq_name).rhof_veff_ratio = rhof_veff_ratio;
    end
    
end

%% Plot output
if parsed_argins.is_plot
    plot_scintillation_realization(out, sim_params.severity, ...
        cspsm_root_dir);
end

%% Play satellite-receiver scenarios
if parsed_argins.is_play
    play(sat_scen);
end

end

