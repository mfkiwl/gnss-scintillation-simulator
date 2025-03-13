function rhof_veff_ratios_L1_struct = get_all_rhof_veff_ratios(training_data_general_params, mask_angle_deg)
% get_all_rhof_veff_ratios
%
% Syntax:
%   rhof_veff_ratios_struct = get_all_rhof_veff_ratios(training_data_general_params, mask_angle_deg)
%
% Description:
%   Computes the mean scaling parameter (rho_F / v_eff) at L1 for every 
%   available satellite (PRN) and a set of simulated ionospheric zonal drift
%   velocities. For each receiver, the function:
%       1) Extracts satellite ephemeris from RINEX files.
%       2) Determines the satellites in line-of-sight (LOS) using a specified 
%          elevation mask.
%       3) Converts the UTC date/time to GPS week and seconds, adjusting for
%          week rollover.
%       4) Generates the receiver trajectory (using GenUserTraj) and computes 
%          satellite geometry parameters via PropGeomCalc.
%       5) Computes the effective IPP range and then applies equations
%          as shown in [1] to obtain the mean ratio (rho_F / v_eff) at L1.
%       6) Organizes the computed ratios into a table for each receiver where
%          the first column ("DriftVelocities") lists the simulated drift 
%          velocities and subsequent columns (named after the PRNs) hold the 
%          corresponding computed parameter (scalar double).
%
% Inputs:
%   training_data_general_params - Struct containing the simulation parameters,
%       including:
%           .date_time           : MATLAB date vector [year, month, day, hour, minute, second]
%           .simulation_time     : Duration of the simulation (seconds)
%           .rx_pos_rad          : Struct with receiver positions in [Latitude, Longitude, Height] (radians)
%           .rx_vel              : Receiver velocity (m/s)
%           .drift_velocities    : Matrix of simulated ionospheric drift velocities (m/s)
%           .ipp_height          : Ionospheric pierce point height (m)
%           .gps_bands           : [L1, L2, L5] frequencies in Hz
%           .c                   : Speed of light (m/s)
%
%   mask_angle_deg - Elevation mask angle in degrees. Satellites with an
%                    elevation below this value are excluded from the LOS 
%                    computation.
%
% Outputs:
%   rhof_veff_ratios_struct - Struct whose fields correspond to receiver
%                             cities. Each field contains a table with:
%                               - Row names defined by the drift velocities.
%                               - A first column "DriftVelocities" listing the drift values.
%                               - Additional columns (named after PRNs) holding the
%                                 computed mean (rho_F / v_eff) at L1 (scalar double).
%
% Notes:
%   - The effective IPP range is computed as: 
%         ipp_range * (sat_receiver_range - ipp_range) / sat_receiver_range,
%     and the mean ratio is computed as:
%         mean( sqrt(effective_ipp_range) / (veff * sqrt((2*pi*f_L1)/c) ) ).
%   - The approach follows that introduced in [1] for scintillation simulation.
%
% Dependencies:
%   UT2GPStime      - Converts UTC date/time to GPS week and seconds.
%   GenUserTraj     - Generates the receiver trajectory (LLH) data.
%   extract_all_rinex_eph - Extracts satellite ephemeris from RINEX files.
%   get_sats_in_los - Determines satellites in LOS based on the elevation mask.
%   PropGeomCalc    - Computes satellite geometry and IPP parameters.
%
% References:
%   [1] Jiao, Yu, Rino, Charles, Morton, Yu (Jade), Carrano, Charles, 
%       "Scintillation Simulation on Equatorial GPS Signals for Dynamic 
%       Platforms," ION GNSS+ 2017, pp. 1644-1657, 
%       https://doi.org/10.33012/2017.15258
%
% See also:
%   UT2GPStime, GenUserTraj, extract_all_rinex_eph, get_sats_in_los, PropGeomCalc
%
% Example:
%   training_data_general_params.date_time = [2023, 01, 10, 12, 00, 00];
%   training_data_general_params.simulation_time = 300;
%   training_data_general_params.rx_pos_rad = struct('CityA', struct('Latitude', rad, 'Longitude', rad, 'Height', m));
%   training_data_general_params.rx_vel = [0, 0, 0];
%   training_data_general_params.drift_velocities = [ ... ]; % e.g., [v_north; v_east; v_vertical]
%   training_data_general_params.ipp_height = 350e3;
%   training_data_general_params.gps_bands = [1.57542e9, 1.22760e9, 1.17645e9];
%   training_data_general_params.c = 299792458;
%   mask_angle_deg = 15;
%   rhof_veff_ratios_struct = get_all_rhof_veff_ratios(training_data_general_params, mask_angle_deg);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    ephs = extract_all_rinex_eph(training_data_general_params.date_time);
    
    % Convert mask angle from degrees to radians.
    mask_angle_rad = mask_angle_deg * pi/180;
    % Determine the satellites in line-of-sight.
    sats_in_los_struct = get_sats_in_los(ephs,training_data_general_params, mask_angle_rad);
    
    % Compute the simulation time series used by PropGeomCalc function (in seconds).
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(training_data_general_params.date_time);
    gps_time_sec_end = gps_time_sec_start + training_data_general_params.simulation_time - 1;
    
    % Build a 2-row array: first row is the GPS week number, second row is seconds.
    gps_week_in_seconds(1, :) = ones(1, training_data_general_params.simulation_time) * gps_week;
    gps_week_in_seconds(2, :) = gps_time_sec_start : gps_time_sec_end;
    
    % Adjust for crossing into a new GPS week (604800 seconds per week).
    index_new_week = find(gps_week_in_seconds(2, :) >= 604800);
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;

    % Initialize the output structure.
    rhof_veff_ratios_L1_struct = struct();

    % Get the names of receiver cities names from `training_data_general_params`.
    cities = fieldnames(training_data_general_params.rx_pos_rad);
    
    % Loop over each receiver city.
    for city_idx = 1:numel(cities)
        % Get the current city field name.
        city_field = cities{city_idx};

         % Build a general parameters struct for this receiver.
        general_parameters = struct();
        general_parameters.rx_pos = [training_data_general_params.rx_pos_rad.(city_field).Latitude, ...
            training_data_general_params.rx_pos_rad.(city_field).Longitude,...
            training_data_general_params.rx_pos_rad.(city_field).Height];
        general_parameters.rx_vel = training_data_general_params.rx_vel;
        general_parameters.simulation_time = training_data_general_params.simulation_time;

        % Compute the fixed receiver position (LLH) from the user trajectory.
        origin_llh = GenUserTraj(general_parameters);

        % Get the PRN names (as a string array) for this receiver's LOS satellites.
        prn_string_array = sats_in_los_struct.(city_field);
        % Extract eastward drift velocities (assuming the second row represents eastward drift).
        eastward_drift_velocities = training_data_general_params.drift_velocities(2, :);

        % Create a numeric table with rows for each drift and columns for each PRN.
        % The table will have an initial column "DriftVelocities" and one column per PRN.
        T_numeric = array2table(zeros(numel(eastward_drift_velocities), numel(prn_string_array)), ...
            'VariableNames', cellstr(prn_string_array));
        rhof_veff_ratio_L1_table = addvars(T_numeric, eastward_drift_velocities(:), 'Before', 1, ...
            'NewVariableNames', 'DriftVelocities');

        % Loop over each satellite PRN.
        for prn_str = prn_string_array
            % Loop over each zonal drift velocity.
            for drift_idx = 1 : size(training_data_general_params.drift_velocities, 2)
                sat_geom = PropGeomCalc( ...
                    gps_week_in_seconds, ...
                    training_data_general_params.date_time, ...
                    ephs.(prn_str), ...
                    origin_llh, ...
                    training_data_general_params.ipp_height, ...
                    training_data_general_params.rx_vel, ...
                    training_data_general_params.drift_velocities(:,drift_idx));
            
                sat_receiver_range = sat_geom.sat_rnge;
                ipp_range = sat_geom.rngp;
                veff = sat_geom.veff;
            
                % Effective IPP range [1, Eq. (13)]
                effective_ipp_range = ipp_range .* (sat_receiver_range - ipp_range) ./ sat_receiver_range;
            
                % Mean ratio (rho_F / v_eff) at L1 [1, Eq. (12)]
                % The usage of the average ratio follows the approach introduced by "Joy" 
                % in the gnss-scintillation-simulator-2-param repository.
                rhof_veff_ratio_L1_table{drift_idx, prn_str} = mean( ...
                    sqrt(effective_ipp_range) ./ (veff * sqrt((2 * pi * training_data_general_params.gps_bands(1)) / training_data_general_params.c)) ...
                );
            end
        end

        % Save the table for this receiver into the output structure.
        rhof_veff_ratios_L1_struct.(city_field) = rhof_veff_ratio_L1_table;
    end
end