function sim_params = set_scenario(log, cpssm_root_dir, ...
    sim_params, parsed_argins)
%GET_SCENARIO Summary of this function goes here
%   Detailed explanation goes here

%% get RINEX file
[sim_params.time_range, rinex] = get_rinex(cpssm_root_dir, parsed_argins);

%% Get line-of-sight (LOS) satellites

% create a satellite scenario
sat_scen = satelliteScenario(sim_params.time_range.start, ...
    sim_params.time_range.end, sim_params.geo_tsamp);

% add all satellites from the RINEX ephemerides in the scenario
satellite(sat_scen, rinex);

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

% add receiver in the scenario as a moving platform
platform(sat_scen, rx_traj, 'Receiver');

% get parameters of the satellites in LOS with the receiver
ac = access(sat_scen.Satellites, sat_scen.Platforms(1));
los_sats_params = ac.accessIntervals;

% set the simulation parameters' frequencies and constellations based on
% the LOS satellites' parameters and user input arguments
sim_params = set_constellation_freq_svid(log, sim_params, ...
    parsed_argins, los_sats_params);

% get user-filtered LOS sats IDS (filtered by constellation and SV IDS)
filtered_los_sats_params = get_filtered_los_sat_params(log, sim_params, ...
    los_sats_params);

% get only the satellites in line-of-sight with the receiver
is_los_sat = ismember(sat_scen.Satellites.Name, filtered_los_sats_params.Source);
for nonlos_sat = sat_scen.Satellites(~is_los_sat)
    delete(nonlos_sat)
end
% set scenario
sim_params.satelliteScenario = sat_scen;
end
