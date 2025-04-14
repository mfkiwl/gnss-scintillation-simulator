function rx_traj_llh = get_rx_traj(general_parameters)
% get_rhof_veff_ratio
%
% Syntax:
%   rx_traj_llh = get_rx_traj(general_parameters)
%
% Description:
%   Computes the receiver tracjectory in LLH (latitude [rad], longitude
%   [rad], height [m]) for all timestamps. It starts the origin LLH for all
%   timestamps. Then, it propagates the receiver trajectory based on
%   receiver velocity.
%
% Inputs:
%   general_parameters - Struct containing the required parameters and sett
%                ings for trajectory and geometry calculations. For this 
%                script, it uses:
%       .simulation_time: Duration for which data is computed (seconds)
%       .rx_pos         : Receiver position in [lat [rad], lon [rad], height [m])
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
% east-west velocity on the earth arc (m/s, eastward +)
eastwest_v = general_parameters.rx_vel(1);
% north-south velocity on the earch arc (m/s, northward +)
northsouth_v = general_parameters.rx_vel(2);
% up-down velocity (m/s, upward +)
updown_v = general_parameters.rx_vel(3);
% simulation time (in seconds)
simulation_time = general_parameters.simulation_time; 
% Earth radius (m)
earthRadius = 6378.137e3;
% receiver tracjectory [LLH x sec]
rx_traj_llh = general_parameters.rx_pos .* ones(3,simulation_time);

%% Main routine
if updown_v ~= 0
    for ii = 2:simulation_time
        rx_traj_llh(3,ii) = rx_traj_llh(3,ii-1)+updown_v;
    end
end
if northsouth_v ~= 0
    northsouth_radv = northsouth_v./(earthRadius+rx_traj_llh(3,:));  % angular velocity (rad/s)
    for ii = 2:simulation_time
        rx_traj_llh(1,ii) = rx_traj_llh(1,ii-1)+northsouth_radv(ii-1);
    end
end
if eastwest_v ~= 0
    eastwest_radv = eastwest_v./(earthRadius+rx_traj_llh(3,:));  % angular velocity (rad/s)
    for ii = 2:simulation_time
        rx_traj_llh(2,ii) = rx_traj_llh(2,ii-1)+eastwest_radv(ii-1);
    end
end