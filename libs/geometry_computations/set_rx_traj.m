function rx = set_rx_traj(rx, origin, t_samp, time_range, earth_radius)
% set_rx_traj
%
% Syntax:
%   sim_params = set_rx_traj(rx, sim_time, earth_radius)
%
% Description:
%   Computes the receiver tracjectory in LLH (latitude [rad], longitude
%   [rad], height [m]) for all timestamps. It begins with origin LLH at
%   the first timestamp. Then, it propagates the receiver trajectory based
%   on receiver (rx) linear and angular velocities. After the computation,
%   the rx trajectory is set at `sim_params.rx` as a new field named
%   `traj`.
%
% Inputs:
%   rx                - Receiver parameters
%   sim_time          - Duration for which data is computed (s)
%   earth_radius      - Earth radius (m)
%   rx_origin        - (optional, [rad, rad, m], 3x1 array)
%                       Receiver position as [latitude; longitude; height].
%                       Default: [0.3876; 1.9942; 59.6780].
%
% Joy
% Written:  03/28/2017
% Modified by Dongyang: 10/18/2018
% Modified by Rodrigo Florindo (https://orcid.org/0000-0003-0412-5583) (06/02/2025):
%   - Substituted the `userInput` from the `gnss-scintillation-simulator-2-
%   param` repository code by the `general_parameters`.
% 
% Modified by Rubem Pacelli (https://orcid.org/0000-0001-5933-8565) (16 Apr 2025):
%   - Heavy refector on the variables. The computation method is only thing
%   kept untouched

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
time_array = time_range.start:seconds(t_samp):time_range.end;
sim_duration = seconds(time_range.end - time_range.start);

%% Propagate recever trajectory - height
traj_height = origin.height + rx.vel.downup * (0:t_samp:sim_duration);

%% Propagate recever trajectory - latitude
% the total radius of the rx trajectory from the center of Earth
rx_traj_radius = earth_radius + traj_height;

% NOTE: in general, if you have a point moving in a circle (or roughly
% circular path) with radius R in metter, the angular velocity w (in
% radians per second) is given by
% w = v_perp / R
% where v_perp is the cross-radial component in m/s (the unique component
% that contributes to the angular velocity)
% SEE: https://en.wikipedia.org/wiki/Angular_velocity#Particle_in_two_dimensions

% angular velocity (rad/s) in the south-north direction
rx_northsouth_angv = rx.vel.southnorth./rx_traj_radius;

% propagate receiver pos in latitude
traj_lat = origin.lat + rx_northsouth_angv .* (0:(t_samp-1));

%% Propagate recever trajectory - logitude

% angular velocity (rad/s) in the west-east direction
rx_westeast_angv = rx.vel.westeast./(rx_traj_radius);
% propagate receiver pos in longitude
traj_long = origin.long + rx_westeast_angv .* (0:(t_samp-1));

%% Set rx trajectory in general parameters
traj = timetable( ...
    time_array.', ...
    traj_height.', ...
    traj_lat.', ...
    traj_long.', ...
    'VariableNames', {'Height','Latitude','Longitude'} ...
);

rx.traj = traj;
