function sim_params = set_sim_params_from_argin(sim_params, parsed_argin)
%% User-defined or dafaulted parameters

%% drift velocity
sim_params.drift_vel_ned = parsed_argin.drift_vel_ned;

%% receiver velocity
% TODO: use the tajectory object as a simulation parameter and use its
% velocity instead the user input
sim_params.rx.vel_ned = parsed_argin.rx_vel_ned;

%% IPP altitude
sim_params.ipp_altitude = parsed_argin.ipp_altitude;

%% Severity
sim_params.severity = parsed_argin.severity;

%% Simulation time
sim_params.sim_time = parsed_argin.sim_time;

%% Sampling time
sim_params.t_samp = parsed_argin.t_samp;
end

