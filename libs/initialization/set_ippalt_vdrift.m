function [sim_params_drift_vel_ned, sim_params_rx, sim_params_ipp_altitude, sim_params_severity] = set_ippalt_vdrift(rx_vel_ned, ...
    drift_vel_ned, ipp_altitude, severity)
%% User-defined or dafaulted parameters

%% drift velocity
sim_params_drift_vel_ned = drift_vel_ned;

%% receiver velocity
% TODO: use the tajectory object as a simulation parameter and use its
% velocity instead the user input
rx.vel_ned = rx_vel_ned;
sim_params_rx = rx;

%% IPP altitude
sim_params_ipp_altitude = ipp_altitude;

%% Severity
sim_params_severity = severity;
end

