function get_rhoF_veff_ratio(general_parameters)
% WRITE THE DOCSTRING HERE
%
% Notes: 
% - I've tried to mantain a similar methodology here from the code
%   `RunPropGeomCalc.m`, in order to get faster results. There are still 
%   many improvements that could achieved here.
%
% Dependencies:
% 

    % Get a time series of the user trajectory
    origin_llh = GenUserTraj(general_parameters);
    % Get the ephemeris of the GPS satellites
    eph = ExtractRINEXeph(general_parameters);

    %% Generate 1-second time samples for the propagation geometry calculation
    
    % UT2GPSTime is a function written by Charles Rino.
    % TODO (Rodrigo): We may substitute this later, since this seems to be very old.
    [gps_time_sec_start,gps_week,~,~]=UT2GPStime(general_parameters.date_time);
    gps_time_sec_end = gps_time_sec_start + general_parameters.simulation_time-1;
    
    % Creates array with (2 x `general_parameters.simulation_time`) elements,
    % where the first line corresponds to the values of the chosen week to
    % perform the simulation and the second line denotes the values in
    % second counting from the start of the week.
    gps_week_in_seconds(1,:) = ones(1,general_parameters.simulation_time)*gps_week;
    gps_week_in_seconds(2,:) = gps_time_sec_start:gps_time_sec_end;
    
    % Handle the case when a week is finishing
    index_new_week = find(gps_week_in_seconds(2,:)>=604800);
    % Corrects the value of the week
    gps_week_in_seconds(1,index_new_week) = gps_week_in_seconds(1,index_new_week)+1;
    % Corrects the value of the time
    gps_week_in_seconds(2,index_new_week) = gps_week_in_seconds(2,index_new_week)-604800;
end