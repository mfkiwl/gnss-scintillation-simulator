function [outputArg1,outputArg2] = cspsm(varargin)
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

%% handle input args
[parsed_input_args, log] = parse_input_args(cspsm_root_dir, varargin{:});

%% get initial simulation parameters
sim_params = get_sim_params(parsed_input_args.rx_vel, ...
    parsed_input_args.drift_vel, parsed_input_args.ipp_alt);

%% get RINEX file

[time_range, rinex] = get_rinex(cspsm_root_dir, ...
    parsed_input_args.is_download_rinex, parsed_input_args.datetime, ...
    parsed_input_args.rinex_filename, parsed_input_args.sim_time);

%% Get line-of-sight (LOS) satellites

% create a satellite scenario
% TODO: At the moment, the sample time is set to one second. As it becomes
% clearer what this value should be, you must reset it as either hardcoded
% or from an input argument
sat_scen = satelliteScenario(time_range.start, time_range.end, 1);

% instantiate satellites for each constellation from the RINEX ephemerides
all_sats = satellite(sat_scen, rinex);

% CAVEAT: It is not clear to me how the velocity is defined or, at least,
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
    [parsed_input_args.rx_origin; parsed_input_args.rx_origin], ...
    [0 seconds(time_range.end - time_range.start)], ...
    'SampleRate', parsed_input_args.t_samp);

rx = platform(sat_scen, rx_traj, 'Receiver');

%play(sat_scen)

% get parameters of the satellites in LOS with the receiver
ac = access(all_sats, rx);

los_stas_params = ac.accessIntervals;

% get user-filtered LOS sats IDS (filtered by constellation and SV IDS)
filtered_los_sats_params = get_filtered_los_sat_params(log, los_stas_params, ...
    parsed_input_args.prn, parsed_input_args.constellation, ...
    time_range);

filtered_los_sats = all_sats(ismember(all_sats.Name, filtered_los_sats_params.Source));

%% Phase screen realization for each satellite-receiver propagation geometry and frequency

% Receiver LLA (latitude, longitude, altitude)
[rx_traj_lla, ~, time_utc] = states(rx, 'CoordinateFrame','geographic');

for filtered_los_sat = filtered_los_sats
    sat_traj_lla = states(filtered_los_sat, 'CoordinateFrame','geographic');
    rhof_veff_ratio = get_rhof_veff_ratio(time_utc, ...
        rx_traj_lla, sim_params.rx.vel, sat_traj_lla, sim_params.drift_vel);
end

end

