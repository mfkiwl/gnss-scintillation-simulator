function out = cspsm(varargin)
% cspsm Compact scintillation phase screen model
%
% Syntax:
%   general_parameters = cspsm()
%   general_parameters = cspsm('rx_origin', value, 'rx_vel', value, ...
%                           'datetime', value, 'prn', value, ...
%                           'sim_time', value, 't_samp', value, ...
%                           'ipp_altitude', value, 'drift_velocity', value)
%
% Description:
%   This function is the entry point of the compact scintillation phase
%   screen model (CSPSM), which models the scintillaiton effects on GNSS
%   receivers.
%
%
%
%
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
% altitude, drift velocity and receiver parameters
% TODO: when adjusting the program for dynamic receiver, set rx parameters
% (trajectory, velocity) directry from the `platform()` instead of the
% parsed input parameters
sim_params = set_sim_params_from_argin(sim_params, parsed_argins);

%% get RINEX file
[sim_params.time_range, rinex] = get_rinex(cspsm_root_dir, parsed_argins);

%% Get line-of-sight (LOS) satellites

% create a satellite scenario
% TODO: At the moment, the sample time is set to one second. As it becomes
% clearer what this value should be, you must reset it as either hardcoded
% or from an input argument
sat_scen = satelliteScenario(sim_params.time_range.start, sim_params.time_range.end, 1);

% instantiate satellites for each constellation from the RINEX ephemerides
all_sats = satellite(sat_scen, rinex);

% FIXME: It is not clear to me how the velocity is defined or, at least,
% maintained constant. `geoTrajectory()` seems to require a starting and
% ending position, as well as the travel duration to compute the
% trajectory. It has some ways to add velocity (or speed) information, but
% it doesn't seem straightforward to set a constant velocity. It is because
% the velocity seems to have a strong dependence on the ending position: A
% distant path seems to make the body accelerate during the trajectory. In
% other words, it seems to have a dynamic velocity. For the sake of first
% try, the receiver is being here to be static
% SEE: https://www.mathworks.com/help/satcom/ref/geotrajectory-system-object.html
rx_traj = geoTrajectory( ...
    [parsed_argins.rx_origin; parsed_argins.rx_origin], ...
    [0 sim_params.sim_time], ...
    'SampleRate', parsed_argins.t_samp);

rx = platform(sat_scen, rx_traj, 'Receiver');

%play(sat_scen)

% get parameters of the satellites in LOS with the receiver
ac = access(all_sats, rx);
los_sats_params = ac.accessIntervals;

% set the simulation parameters' frequencies and constellations based on
% the LOS satellites' parameters and user input arguments
sim_params = set_constellation_freq_svid(log, sim_params, ...
    parsed_argins, los_sats_params);

% get user-filtered LOS sats IDS (filtered by constellation and SV IDS)
filtered_los_sats_params = get_filtered_los_sat_params(log, sim_params, los_sats_params);

filtered_los_sats = all_sats(ismember(all_sats.Name, ...
    filtered_los_sats_params.Source));

%% Phase screen realization for each satellite-receiver propagation geometry and frequency

% get receiver LLA (latitude, longitude, altitude) and UTC time
[rx_traj_lla, ~, time_utc] = states(rx, 'CoordinateFrame','geographic');

% TODO: get the receiver velocity from `states(rx, ...)` instead of the
% user inputs
rx_vel_ned = repmat(sim_params.rx.vel_ned.', 1, numel(time_utc));

seed = 0;

% initialize output
for constellation = sim_params.constellations
    out.(constellation) = struct([]); % 0x0 struct
end

% for all filtered LOS sat (of a given contellation)
for i = 1:numel(filtered_los_sats)
    % get LOS sat and constellation name
    filtered_los_sat = filtered_los_sats(i);
    sat_constellation = filtered_los_sat.OrbitPropagator;
    % save sat and receiver to output
    out.(sat_constellation)(end+1).sat = filtered_los_sat;
    out.(sat_constellation)(end).rx = rx;
    % set satellite trajectory and velocities
    [sat_traj_lla, sat_vel_ned, ~]= states(filtered_los_sat, 'CoordinateFrame','geographic');
    % for this geometry propagation, get the ρF/veff for the reference
    % frequency
    % TODO: shrink argins adding `rx` and `filtered_los_sat` within
    rhof_veff_ratio_ref = get_rhof_veff_ratio(time_utc, ...
        rx_traj_lla, rx_vel_ned, sat_traj_lla, sat_vel_ned, ...
        sim_params.drift_vel_ned, sim_params.ipp_altitude, ...
        sim_params.cte.spectral.freq_ref.value, sim_params.cte.c);

    % for all valid frequencies for this contellation
    for j = 1:numel(sim_params.freqs.(sat_constellation).name)
        freq_name = sim_params.freqs.(sat_constellation).name(j);
        freq_value = sim_params.freqs.(sat_constellation).value(j);

        % extrapolate ρF/veff, U, and μ₀ from `freq_ref` to `freq`
        [sim_params.rhof_veff_ratio, sim_params.spectral] = ...
            freq_extrapolate(sim_params.cte.spectral.(sim_params.severity), ...
            rhof_veff_ratio_ref, sim_params.cte.spectral.freq_ref.value, freq_value);
        
        % TODO: shrink argins
        [scint_field, norm_phase_psd, detrended_phase, mu, doppler_frequency] = ...
            get_scintillation_time_series(sim_params.sim_time, ...
            sim_params.t_samp, ...
            sim_params.spectral, ...
            sim_params.rhof_veff_ratio, ...
            seed);

        % compute amplitude and phase time series of the received
        % scintillation signal
        amplitude = abs(scint_field.');
        phase = unwrap(angle(scint_field.'));
        
        % outputs
        out.(sat_constellation)(end).(freq_name).scint_field = scint_field;
        out.(sat_constellation)(end).(freq_name).amplitude = amplitude;
        out.(sat_constellation)(end).(freq_name).phase = phase;
        out.(sat_constellation)(end).(freq_name).detrended_phase = detrended_phase;
        out.(sat_constellation)(end).(freq_name).mu = mu;
        out.(sat_constellation)(end).(freq_name).doppler_frequency = doppler_frequency;
        out.(sat_constellation)(end).(freq_name).norm_phase_psd = norm_phase_psd;
        out.(sat_constellation)(end).(freq_name).spectral_params = sim_params.spectral;
        out.(sat_constellation)(end).(freq_name).rhof_veff_ratio = sim_params.rhof_veff_ratio;
    end
    
end

%% Plot output

if parsed_argins.is_plot
    plot_scintillation_realization(out, sim_params.severity, cspsm_root_dir);
end

end

