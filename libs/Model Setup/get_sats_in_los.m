function sats_in_los_struct = get_sats_in_los(ephs, training_data_general_params, mask_angle)
% get_sats_in_los
%
% Syntax:
%   sats_in_los_struct = get_sats_in_los(ephs, training_data_general_params, mask_angle)
%
% Description:
%   Determines which satellites are in line-of-sight (LOS) for each receiver.
%   For each receiver position, the function computes satellite positions at
%   the simulation start and end times (accounting for GPS week rollover),
%   converts them to a topocentric coordinate system relative to the receiver,
%   and checks if their elevation is above a specified mask angle.
%
% Inputs:
%   ephs - Struct containing satellite ephemeris data. Each field should be
%          named by a satellite PRN.
%   training_data_general_params - Struct with simulation parameters, including:
%          .date_time   : [YYYY MM DD hh mm ss] vector of simulation start.
%          .simulation_time : Total simulation duration in seconds.
%          .rx_pos_rad  : Struct with receiver positions in radians. Each
%                         field corresponds to a receiver and contains
%                         subfields 'Latitude', 'Longitude', and 'Height'.
%   mask_angle - Minimum elevation angle (in radians) required for a satellite
%                to be considered in LOS.
%
% Outputs:
%   sats_in_los_struct - Struct where each field corresponds to a receiver.
%                        Each receiver field contains a string array of satellite
%                        PRNs that are in LOS (above the mask angle).
%
% Example:
%   sats_struct = get_sats_in_los(ephs, training_data_general_params, deg2rad(15));
%
% Dependencies:
%   UT2GPStime, satposvel, ecf2llhT, llh2tcsT.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Compute the simulation start and end GPS times (in seconds).
    [gps_time_sec_start, ~, ~, ~] = UT2GPStime(training_data_general_params.date_time);
    gps_time_sec_end = gps_time_sec_start + training_data_general_params.simulation_time - 1;
    
    gps_week_in_seconds = [gps_time_sec_start, gps_time_sec_end];
    
    % Correct for crossing into a new GPS week (604800 seconds).
    week_in_seconds = 604800;
    idx_new_week = gps_week_in_seconds >= week_in_seconds;
    gps_week_in_seconds(idx_new_week) = gps_week_in_seconds(idx_new_week) - week_in_seconds;
    
    % Get the names of receiver cities (as cell array of char).
    cities = fieldnames(training_data_general_params.rx_pos_rad);
    
    % Get the list of satellite PRNs (as cell array of char).
    prn_list = fieldnames(ephs);
    
    % Initialize the output struct.
    sats_in_los_struct = struct();
    
    % Loop over each receiver (city).
    for city_idx = 1:numel(cities)
        % Get the current city field name.
        city_field = cities{city_idx};
        
        % Extract the receiver's position (in radians) from the general params.
        city_pos = training_data_general_params.rx_pos_rad.(city_field);
        origin_llh = [city_pos.Latitude; city_pos.Longitude; city_pos.Height];
        
        % Preallocate available_prns to maximum possible size (number of PRNs)
        available_prns = strings(numel(prn_list), 1);
        count = 0;
        
        % Loop over each satellite PRN.
        for prn_idx = 1:numel(prn_list)
            prn_field = prn_list{prn_idx};
            
            % Obtain the satellite position (in ECF) at the simulation start
            % and end times.
            [start_sat_pos_ecf, ~] = satposvel(gps_week_in_seconds(1), ephs.(prn_field));
            [end_sat_pos_ecf, ~]   = satposvel(gps_week_in_seconds(2), ephs.(prn_field));
            
            % Convert satellite positions from ECF to geodetic (LLH) coordinates.
            start_sat_llh = ecf2llhT(start_sat_pos_ecf);
            end_sat_llh   = ecf2llhT(end_sat_pos_ecf);
            
            % Convert LLH coordinates to Topocentric Coordinate System (TCS)
            % at the receiver's location.
            start_sat_tcs = llh2tcsT(start_sat_llh, origin_llh);
            end_sat_tcs   = llh2tcsT(end_sat_llh, origin_llh);
            
            % Compute elevation angles using the TCS coordinates.
            start_sat_elev = atan2(start_sat_tcs(3), sqrt(start_sat_tcs(1).^2 + start_sat_tcs(2).^2));
            end_sat_elev   = atan2(end_sat_tcs(3), sqrt(end_sat_tcs(1).^2 + end_sat_tcs(2).^2));
            
            % Check if the satellite is above the mask angle at both times.
            if ((start_sat_elev > mask_angle) && (end_sat_elev > mask_angle))
                count = count + 1;
                available_prns(count) = string(prn_field);
            end
        end
        
        % Trim the available_prns array to the actual count.
        available_prns = available_prns(1:count);
        
        % Store the list of available PRNs for this receiver.
        sats_in_los_struct.(city_field) = available_prns;
    end
end
