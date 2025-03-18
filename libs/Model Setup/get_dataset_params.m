function data_set_params = get_dataset_params(drift_velocities_amount, carrier_to_noise_ratios_amount, cities_amount, mask_angle_deg)
% get_dataset_params
%
% Syntax:
%   data_set_params = get_dataset_params(drift_velocities_amount, carrier_to_noise_ratios_amount, mask_angle_deg)
%
% Description:
%   Generates a struct with simulation parameters organized into two top‐level
%   fields: 'general' and 'specific'. 
%
%   general includes:
%       c, gps_bands, date_time, dt, simulation_time, ipp_height, drift_velocities,
%       carrier_to_noise_ratios, and irr_params.
%
%   specific includes, for each city:
%       rx_init_llh, rx_velocity, prn_string_array, and rhof_veff_table.
%
% Dependencies:
%   get_irregularity_parameters, local_extrapolate, dms_str_2_rad,
%   extract_all_rinex_eph, UT2GPStime, get_sats_in_los, GenUserTraj,
%   PropGeomCalc.
%
% Inputs:
%   drift_velocities_amount        - Number of drift velocity values to generate.
%   carrier_to_noise_ratios_amount - Number of carrier-to-noise ratio values.
%   cities_amount                  - Amount of cities to be considered
%   mask_angle_deg                 - Elevation mask angle in degrees.
%
% Outputs:
%   data_set_params - Struct with the fields 'general' and 'specific'.
%
% Notes:
%   - DMS strings use the format: deg°min''sec"Direction (N,S,E,W).
%   - Function dms_str_2_rad must be available in the MATLAB path.
%
% References:
%  [1] https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/f46820626881276a4bfb6d4dcf1225f74e6cd372/Libraries/PropGeomCalc/RunPropGeomCalc.m#L48C1-L48C13
%  [2] https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/f46820626881276a4bfb6d4dcf1225f74e6cd372/Libraries/PropGeomCalc/RunPropGeomCalc.m#L49
%  [3] Y. Jiao, D. Xu, C. L. Rino, Y. T. Morton and C. S. Carrano, "A 
%      Multifrequency GPS Signal Strong Equatorial Ionospheric 
%      Scintillation Simulator: Algorithm, Performance, and 
%      Characterization," in IEEE Transactions on Aerospace and Electronic 
%      Systems, vol. 54, no. 4, pp. 1947-1965, Aug. 2018, 
%      doi: 10.1109/TAES.2018.2805232
% [4] Carrano, C. S., and C. L. Rino (2016), A theory of scintillation for 
%     two-component power law irregularity spectra: Overview and numerical
%     results, Radio Sci., 51, 789–813, doi:10.1002/2015RS005903.
% [5] https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/f46820626881276a4bfb6d4dcf1225f74e6cd372/Libraries/GenScintFieldRealization/RunGenScintFieldRealization.m#L32
% [6] https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/f46820626881276a4bfb6d4dcf1225f74e6cd372/Libraries/GenScintFieldRealization/RunGenScintFieldRealization.m#L49
% [7] https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/f46820626881276a4bfb6d4dcf1225f74e6cd372/Libraries/GenScintFieldRealization/RunGenScintFieldRealization.m#L59
% [8] Jiao, Y., Rino, C., and Morton, Y. T. (2018) Ionospheric 
%     Scintillation Simulation on Equatorial GPS Signals for Dynamic 
%     Platforms. Journal of the Institute of Navigation, 65: 263–274. 
%     doi: https://doi.org/10.1002/navi.231.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    %% General Parameters
    general = struct();
    
    % Simulation constants
    general.c = 299792458;  % Speed of light in vacuum (m/s)
    
    % GPS frequencies as multiples of the 10.23 MHz fundamental frequency.
    general.gps_bands = [154*10.23e6, 120*10.23e6, 115*10.23e6];  % L1, L2, L5 in Hz
    
    % Double array with the start date in the format [YYYY MM DD hh mm ss].
    general.date_time = [2014 01 02 10 00 00];  
    % Sampling interval (s).
    general.dt = 0.01;
    general.simulation_time = 300;
    % Ionospheric pierce point height (m). NOTE: This is the same as used
    % in the original version of the compact GNSS scintillation simulator [1].
    general.ipp_height = 350000;
    
    % Linearly spaced values from 25 to 125 m/s, with the amount of samples
    % being configured by `drift_velocities_amount`. NOTE: The bounds
    % [25,125] m/s for the zonal drift velocity comes from [3, Section IV.D]
    drift_velocities = linspace(25,125,drift_velocities_amount);

    % Input the drift velocities in a (3 × `drift_velocities_amount`) matrix. 
    % The first and third row corresponds to the downward and southward
    % direction, respectively. NOTE: We assume that only a eastward zonal
    % drift velocity happens, as it is done in [2].
    general.drift_velocities = [zeros(1,drift_velocities_amount); drift_velocities; zeros(1,drift_velocities_amount)];
    general.drift_velocities_amount = drift_velocities_amount;
    
    % Carrier-to-noise ratios linearly spaced from 25 to 45 dB-Hz.
    general.carrier_to_noise_ratios = linspace(25,45,carrier_to_noise_ratios_amount);
    general.carrier_to_noise_ratios_amount = carrier_to_noise_ratios_amount;

    % Here, a beforehand configured set of the irregularity parameters 
    % I = {U, mu0, p1, p2} is obtained for the L1 frequency band. NOTE: 
    % The irregularity parameter values for the severity classes
    % {weak, moderate and strong} were identified based on [4].
    L1_irr_params_set = get_irregularity_parameters();

    % Here, the function `local_extrapolate` is used to extrapolate the
    % values of the irregularity parameters for L1 to L2 and L5,
    % respectively.
    severity_names = ["strong", "moderate", "weak"];
    for severity_name = severity_names
        [U_L2, mu0_L2] = local_extrapolate(L1_irr_params_set.(severity_name).U, ...
            L1_irr_params_set.(severity_name).mu0, ...
            L1_irr_params_set.(severity_name).p1, ...
            L1_irr_params_set.(severity_name).p2, ...
            general.gps_bands(1), ...
            general.gps_bands(2));
        [U_L5, mu0_L5] = local_extrapolate(L1_irr_params_set.(severity_name).U, ...
            L1_irr_params_set.(severity_name).mu0, ...
            L1_irr_params_set.(severity_name).p1, ...
            L1_irr_params_set.(severity_name).p2, ...
            general.gps_bands(1), ...
            general.gps_bands(3));
        general.irr_params.(severity_name).p1 = L1_irr_params_set.(severity_name).p1;
        general.irr_params.(severity_name).p2 = L1_irr_params_set.(severity_name).p2;
        general.irr_params.(severity_name).U_L1 = L1_irr_params_set.(severity_name).U;
        general.irr_params.(severity_name).mu0_L1 = L1_irr_params_set.(severity_name).mu0;
        general.irr_params.(severity_name).U_L2 = U_L2;
        general.irr_params.(severity_name).mu0_L2 = mu0_L2;
        general.irr_params.(severity_name).U_L5 = U_L5;
        general.irr_params.(severity_name).mu0_L5 = mu0_L5;
    end

    %% Specific Parameters
    % Define receiver positions (in DMS) for each city and their velocities.
    % NOTE: 
    % 1) For now, i've arbitrarily selected 13 different places on earth
    %    so we can obtain a meaningful set of values of the time scaling 
    %    parameter based on realistic satellite-receiver range dynamics.
    % 2) I've used Google Earth Pro to obtain the LLH of each city. 
    %    It showed the coordinates the format D°mm's.ss"(NSWE), which i
    %    needed to refactor it into two separate strings 
    %    as 'D°mm'' 's.ss"(NSWE)', for correct formatting here.
    receivers_positions_dms = struct( ...
        'ascension_island',    struct("Latitude", '7°57''4.97"S',   "Longitude", '14°21''1.06"W',   "Height", 728, 'Velocity', [0;0;0]), ...
        'hong_kong',           struct("Latitude", '22°19''2.21"N',  "Longitude", '114°10''2.33"E',  "Height", 50, 'Velocity', [0;0;0]), ...
        'fortaleza',           struct("Latitude", '3°44''22.26"S',  "Longitude", '38°31''21.46"W',  "Height", 26, 'Velocity', [0;0;0]), ...
        'chiang_mai',          struct("Latitude", '18°47''6.09"N',  "Longitude", '98°59''15.70"E',  "Height", 310, 'Velocity', [0;0;0]), ...
        'port_moresby',        struct("Latitude", '9°28''47.83"S',  "Longitude", '147°9''0.58"E',   "Height", 24, 'Velocity', [0;0;0]), ...
        'san_antonio',         struct("Latitude", '12°56''39.09"N', "Longitude", '121°29''30.07"E', "Height", 12, 'Velocity', [0;0;0]), ...
        'pontianak',           struct("Latitude", '0°1''45.05"S',   "Longitude", '109°20''21.23"E', "Height", 3, 'Velocity', [0;0;0]), ...
        'thiruvananthapuram',  struct("Latitude", '8°31''12.61"N',  "Longitude", '76°55''58.58"E',  "Height", 55, 'Velocity', [0;0;0]), ...
        'nairobi',             struct("Latitude", '1°17''45.22"S',  "Longitude", '36°49''8.98"E',   "Height", 1680, 'Velocity', [0;0;0]), ...
        'abuja',               struct("Latitude", '9°3''4.19"N',    "Longitude", '7°30''1.37"E',    "Height", 498, 'Velocity', [0;0;0]), ...
        'quito',               struct("Latitude", '0°13''23.50"S',  "Longitude", '78°30''48.39"W',  "Height", 2942, 'Velocity', [0;0;0]), ...
        'iquitos',             struct("Latitude", '3°44''37.55"S',  "Longitude", '73°15''5.59"W',   "Height", 95, 'Velocity', [0;0;0]), ...
        'sao_jose_dos_campos', struct("Latitude", '23°13''11.37"S', "Longitude", '45°53''29.72"W',  "Height", 603, 'Velocity', [0;0;0]) ...
    );

    % Check that the requested number of cities does not exceed the available count.
    available_cities = fieldnames(receivers_positions_dms);
    if cities_amount > numel(available_cities)
        error('Only %d cities are available, but cities_amount = %d was requested.', numel(available_cities), cities_amount);
    end
    % Use only the first 'cities_amount' cities.
    selected_cities = available_cities(1:cities_amount);

    % Convert the DMS positions to radians.
    % NOTE: the function dms_str_2_rad is a function that i coded by myself.
    receivers_positions_rad = structfun(@(p) struct( ...
        'Latitude',  dms_str_2_rad(p.Latitude), ...
        'Longitude', dms_str_2_rad(p.Longitude), ...
        'Height',    p.Height, ...
        'Velocity',  p.Velocity), receivers_positions_dms, 'UniformOutput', false);
    
    % Extract all GPS ephemeris available in the RINEX navigation file.
    ephs = extract_all_rinex_eph(general.date_time);
    
    % Compute the GPS time series (week & seconds) used for the `PropGeomCalc` 
    % function developed by the Colorado Boulder University group.
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(general.date_time);
    gps_time_sec_end = gps_time_sec_start + general.simulation_time - 1;
    
    gps_week_in_seconds = zeros(2, general.simulation_time);
    % The first row corresponds to the week of the simulation
    gps_week_in_seconds(1, :) = ones(1, general.simulation_time) * gps_week;
    % The second row refers to the exact time in seconds of the simulation,
    % computed in reference to the beginning of the current week.
    gps_week_in_seconds(2, :) = gps_time_sec_start : gps_time_sec_end;

    % Handles the week changing case.
    index_new_week = find(gps_week_in_seconds(2, :) >= 604800);
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;
    
    specific = struct();
    
    for i = 1:length(selected_cities)
        city = selected_cities{i};

        % Receiver initial position and velocity for this city.
        rx_init_llh = [receivers_positions_rad.(city).Latitude; ...
                       receivers_positions_rad.(city).Longitude; ...
                       receivers_positions_rad.(city).Height];
        rx_velocity = receivers_positions_rad.(city).Velocity;
        
        % Get the in-line-of-sight GPS satellites (PRNs) for this receiver.
        prn_string_array = get_sats_in_los(ephs, rx_init_llh, rx_velocity, general.date_time, general.simulation_time, mask_angle_deg);
        
        % Build a struct for propagating the user position using the 
        % `GenUserTraj` function developed by the Colorado Boulder group.
        % NOTE: This function computes the time series of the user position
        % in latitude, longitude and height (LLH), based on a initial
        % position and its instantaneous velocity.
        user_traj_data = struct();
        user_traj_data.rx_pos = rx_init_llh.';  % ensure row vector
        user_traj_data.rx_vel = rx_velocity;
        user_traj_data.simulation_time = general.simulation_time;
        
        % Generate the receiver trajectory.
        origin_llh = GenUserTraj(user_traj_data);
        
        % Extract eastward drift velocities (assumed to be in row 2).
        eastward_drift_velocities = general.drift_velocities(2, :);
        
       % Create a cell array filled with empty structs.
        empty_struct_cell = repmat( ...
            {struct('rhof_veff_ratio_L1', [], 'rhof_veff_ratio_L2', [], 'rhof_veff_ratio_L5', [])}, ...
            numel(eastward_drift_velocities), ...
            numel(prn_string_array));
        
        % Convert the cell array to a table with column names from prn_string_array.
        T_cell = cell2table(empty_struct_cell, 'VariableNames', cellstr(prn_string_array));
        
        % Add the DriftVelocities column to the table.
        rhof_veff_ratio_table = addvars(T_cell, eastward_drift_velocities(:), 'Before', 1, ...
            'NewVariableNames', 'DriftVelocities');
        
        % Loop over each satellite PRN in the LOS.
        for prn_str = prn_string_array
            % For each simulated drift velocity, compute the geometry parameters.
            for drift_idx = 1:size(general.drift_velocities, 2)
                sat_geom = PropGeomCalc( ...
                    gps_week_in_seconds, ...
                    general.date_time, ...
                    ephs.(prn_str), ...
                    origin_llh, ...
                    general.ipp_height, ...
                    rx_velocity, ...  % use city-specific receiver velocity
                    general.drift_velocities(:, drift_idx));
                
                sat_receiver_range = sat_geom.sat_rnge;
                ipp_range = sat_geom.rngp;
                veff = sat_geom.veff;
                
                % Compute effective IPP range.
                effective_ipp_range = ipp_range .* (sat_receiver_range - ipp_range) ./ sat_receiver_range;
                
                % Compute the mean ratio (rho_F / v_eff) at L1. The method
                % for calculating the time scaling parameter
                % `rhof_veff_ratio_L1` is the same proposed in [5,6,7].
                rhof_veff_ratio_L1 = mean( sqrt(effective_ipp_range) ./ ( veff * sqrt( (2 * pi * general.gps_bands(1) ) / general.c) ) );

                % Scale to obtain the ratios for L2 and L5. 
                % based on [8, Equation 18].
                rhof_veff_ratio_L2 = rhof_veff_ratio_L1 * sqrt(general.gps_bands(1) / general.gps_bands(2));
                rhof_veff_ratio_L5 = rhof_veff_ratio_L1 * sqrt(general.gps_bands(1) / general.gps_bands(3));

                rhof_veff_ratio_table{drift_idx, prn_str}.rhof_veff_ratio_L1 = rhof_veff_ratio_L1;
                rhof_veff_ratio_table{drift_idx, prn_str}.rhof_veff_ratio_L2 = rhof_veff_ratio_L2;
                rhof_veff_ratio_table{drift_idx, prn_str}.rhof_veff_ratio_L5 = rhof_veff_ratio_L5;
            end
        end
        
        % Insert the computed values into the specific parameters for this city.
        specific.(city) = struct(...
            'rx_init_llh', rx_init_llh, ...
            'rx_velocity', rx_velocity, ...
            'prn_string_array', prn_string_array, ...
            'rhof_veff_table', rhof_veff_ratio_table);
    end

    % Combine General and Specific Parameters into the final struct that
    % the function output
    data_set_params = struct();
    data_set_params.general = general;
    data_set_params.specific = specific;
end
