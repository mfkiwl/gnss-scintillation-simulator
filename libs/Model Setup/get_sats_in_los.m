function prn_string_array = get_sats_in_los(ephs, receiver_init_llh, receiver_velocity, datetime, sim_time, mask_angle_deg)

    % Convert mask angle from degrees to radians.
    mask_angle_rad = mask_angle_deg * pi/180;
    
    % Compute the simulation start and end GPS times (in seconds).
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(datetime);
    gps_time_sec_end = gps_time_sec_start + sim_time - 1;
    
    % Build a 2-row array: first row is the GPS week number, second row is seconds.
    gps_week_in_seconds(1, :) = ones(1, 2) * gps_week;
    gps_week_in_seconds(2, :) = [gps_time_sec_start, gps_time_sec_end];
    
    % Correct for crossing into a new GPS week (604800 seconds).
    week_in_seconds = 604800;

    % Adjust for crossing into a new GPS week (604800 seconds per week).
    index_new_week = find(gps_week_in_seconds(2, :) >= week_in_seconds);
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;
    
    % Generate the receiver position evolution across the simulation time.
    user_traj_data = struct();
    user_traj_data.rx_pos = receiver_init_llh;
    user_traj_data.rx_vel = receiver_velocity;
    user_traj_data.sim_time = sim_time;
    rx_traj_llh = set_rx_traj(user_traj_data);

    % Get the list of satellite PRNs (as a cell array of char).
    prn_list = fieldnames(ephs);
    
    % Preallocate available_prns to maximum possible size.
    available_prns = strings(1, numel(prn_list));
    count = 0;
    
    % Loop over each satellite PRN.
    for prn_idx = 1:numel(prn_list)
        prn_field = prn_list{prn_idx};
        
        % Obtain the satellite position (in ECF) at the simulation start and end times.
        [start_sat_pos_ecf, ~] = satposvel(gps_week_in_seconds(2,1), ephs.(prn_field));
        [end_sat_pos_ecf, ~]   = satposvel(gps_week_in_seconds(2,2), ephs.(prn_field));
        
        % Convert satellite positions from ECF to geodetic (LLH) coordinates.
        start_sat_llh = ecf2llhT(start_sat_pos_ecf);
        end_sat_llh   = ecf2llhT(end_sat_pos_ecf);
        
        % Convert LLH coordinates to Topocentric Coordinate System (TCS) relative to the receiver.
        start_sat_tcs = llh2tcsT(start_sat_llh, rx_traj_llh(:,1));
        end_sat_tcs   = llh2tcsT(end_sat_llh, rx_traj_llh(:,2));
        
        % Compute elevation angles using the TCS coordinates.
        start_sat_elev = atan2(start_sat_tcs(3), sqrt(start_sat_tcs(1).^2 + start_sat_tcs(2).^2));
        end_sat_elev   = atan2(end_sat_tcs(3), sqrt(end_sat_tcs(1).^2 + end_sat_tcs(2).^2));
        
        % Check if the satellite is above the mask angle at both times.
        if (start_sat_elev > mask_angle_rad) && (end_sat_elev > mask_angle_rad)
            count = count + 1;
            available_prns(count) = string(prn_field);
        end
    end
    
    % Trim the available_prns array to the actual count.
    prn_string_array = available_prns(1:count);
end
