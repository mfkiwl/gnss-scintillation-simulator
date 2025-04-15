function rx_traj_llh = get_rx_traj(general_parameters)
% get_rhof_veff_ratio
%
% Syntax:
%   rx_traj_llh = get_rx_traj(general_parameters)
%
% Description:
%   Computes the receiver tracjectory in LLH (latitude [rad], longitude
%   [rad], height [m]) for all timestamps. It begins with starting LLH at
%   the first timestamp. Then, it propagates the receiver trajectory based
%   on receiver linear and angular velocities.
%
% Inputs:
%   general_parameters - Struct containing the required parameters and sett
%                ings for trajectory and geometry calculations. For this 
%                script, it uses:
%       .simulation_time: Duration for which data is computed (seconds)
%       .rx_pos         : Receiver position in [lat [rad], long [rad], height [m]]
%       .rx_vel         : Receiver velocity (m/s)
% The components and units for rx_traj_llh are latitude (rad), longitude
% (rad), and height (meter).
% The velocity vector consists of components in east-west, north-south, up-down
% directions.
% Joy
% Written:  03/28/2017
% Modified by Dongyang: 10/18/2018
% Modified by Rodrigo Florindo (06/02/2025):
%   - Substituted the `userInput` from the `gnss-scintillation-simulator-2-
% param` repository code by the `general_parameters`.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
% rx east-west velocity on the earth arc (m/s, eastward +)
rx_eastwest_v = general_parameters.rx_vel(1);
% rx north-south velocity on the earch arc (m/s, northward +)
rx_northsouth_v = general_parameters.rx_vel(2);
% rx up-down velocity (m/s, upward +)
rx_updown_v = general_parameters.rx_vel(3);
% rx latitude
rx_lat = general_parameters.rx_pos(1);
% rx longitude
rx_long = general_parameters.rx_pos(2);
% rx height
rx_height = general_parameters.rx_pos(3);
% simulation time (in seconds)
simulation_time = general_parameters.simulation_time; 
% Earth radius (m)
earth_radius = 6378.137e3;
% receiver tracjectory [LLH x sec]
rx_traj_llh = zeros(3,simulation_time);

%% Propagate recever trajectory - height
rx_traj_llh(3,:) = rx_height + rx_updown_v * (0:(simulation_time-1));

%% Propagate recever trajectory - latitude
% the total radius of the rx trajectory from the center of Earth
rx_total_radius = earth_radius + rx_traj_llh(3,:);

% NOTE: in general, if you have a point moving in a circle (or roughly
% circular path) with radius R in metter, the angular velocity w (in
% radians per second) is given by
% w = v_perp / R
% where v_perp is the cross-radial component in m/s (the unique component
% that contributes to the angular velocity)
% SEE: https://en.wikipedia.org/wiki/Angular_velocity#Particle_in_two_dimensions

% angular velocity (rad/s) in the south-north direction
rx_northsouth_angv = rx_northsouth_v./rx_total_radius;

% propagate receiver pos in latitude
rx_traj_llh(1,:) = rx_lat + rx_northsouth_angv .* (0:(simulation_time-1));

%% %% Propagate recever trajectory - logitude

% angular velocity (rad/s) in the east-west direction
rx_eastwest_angv = rx_eastwest_v./(rx_total_radius);
% propagate receiver pos in longitude
rx_traj_llh(2,:) = rx_long + rx_eastwest_angv .* (0:(simulation_time-1));
