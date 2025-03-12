function training_data_geom_params = get_training_data_general_params(drift_velocities_amount)
% get_training_data_general_params
%
% Syntax:
%   params = get_training_data_geom_params(drift_velocities_amount)
%
% Description:
%   Generates a struct with simulation parameters and receiver geometry.
%   The output struct contains constants, GPS frequencies, receiver
%   positions in both DMS and radian formats, receiver velocity, satellite
%   PRN list, simulation time settings, and ionospheric drift velocities.
%
% Inputs:
%   drift_velocities_amount - Number of drift velocity values to generate.
%
% Outputs:
%   training_data_geom_params - Struct with simulation parameters:
%       .c               : Speed of light in vacuum (m/s).
%       .gps_bands       : GPS frequencies in Hz (L1, L2, L5). The values are
%                          multiples of the 10.23 MHz base frequency.
%       .rx_pos_dms      : Receiver positions in DMS format. Latitude and
%                          Longitude are DMS strings; Height in meters.
%       .rx_pos_rad      : Receiver positions in radians. Latitude and
%                          Longitude in radians; Height in meters.
%       .rx_vel          : Receiver velocity vector [0; 0; 0] (m/s). The
%                          station is assumed to be stationary.
%       .date_time       : Simulation start time as [YYYY MM DD hh mm ss].
%       .simulation_time : Total simulation duration in seconds.
%       .dt              : Sampling interval in seconds.
%       .ipp_height      : Ionospheric pierce point height in meters.
%       .drift_velocity  : 3xn array of drift velocities. The second row has
%                          values linearly spaced from 25 to 125 m/s, with the
%                          first and third rows equal to zero.
%
% Notes:
%   - DMS strings use the format: deg°min''sec"Direction (N,S,E,W).
%   - Function dms_str_2_rad must be available in the MATLAB path.
%
% Example:
%   params = get_training_data_geom_params(10);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    training_data_geom_params = struct();

    % Constant simulation parameters:
    training_data_geom_params.c = 299792458;  % Speed of light in vacuum (m/s)
    
    % GPS frequencies: defined as multiples of the 10.23 MHz fundamental frequency.
    training_data_geom_params.gps_bands = [154*10.23e6, ...  % L1 frequency in Hz
                                    120*10.23e6, ...  % L2 frequency in Hz
                                    115*10.23e6];     % L5 frequency in Hz

    % struct with the receiver_positions in latitude, longitude and height
    % (LLH). Latitude and Longitude is expressed here in degree, minute and
    % seconds (DMS), and height is in meters.
    receivers_positions_dms = struct( ...
        'ascension_island',    struct("Latitude", '7°57''4.97"S',   "Longitude", '14°21''1.06"W',   "Height", 728), ...
        'hong_kong',           struct("Latitude", '22°19''2.21"N',  "Longitude", '114°10''2.33"E',  "Height", 50), ...
        'fortaleza',           struct("Latitude", '3°44''22.26"S',  "Longitude", '38°31''21.46"W',  "Height", 26), ...
        'chiang_mai',          struct("Latitude", '18°47''6.09"N',  "Longitude", '98°59''15.70"E',  "Height", 310), ...
        'port_moresby',        struct("Latitude", '9°28''47.83"S',  "Longitude", '147°9''0.58"E',   "Height", 24), ...
        'san_antonio',         struct("Latitude", '12°56''39.09"N', "Longitude", '121°29''30.07"E', "Height", 12), ...
        'pontianak',           struct("Latitude", '0°1''45.05"S',   "Longitude", '109°20''21.23"E', "Height", 3), ...
        'thiruvananthapuram',  struct("Latitude", '8°31''12.61"N',  "Longitude", '76°55''58.58"E',  "Height", 55), ...
        'nairobi',             struct("Latitude", '1°17''45.22"S',  "Longitude", '36°49''8.98"E',   "Height", 1680), ...
        'abuja',               struct("Latitude", '9°3''4.19"N',    "Longitude", '7°30''1.37"E',    "Height", 498), ...
        'quito',               struct("Latitude", '0°13''23.50"S',  "Longitude", '78°30''48.39"W',  "Height", 2942), ...
        'iquitos',             struct("Latitude", '3°44''37.55"S',  "Longitude", '73°15''5.59"W',   "Height", 95), ...
        'sao_jose_dos_campos', struct("Latitude", '23°13''11.37"S', "Longitude", '45°53''29.72"W',  "Height", 603) ...
    );
    
    % Construct another struct (`receivers_positions_rad`) where latitude and 
    % longitude is expressed in radians instead of DMS. The approach used here
    % iterates through each field of `receiver_positions_dms` and constructs
    % another struct by converting their DMS values to radians and repeating
    % the value of height in meters. The flag 'UniformOutput' is configured to
    % false to ensure that the outputed variable is in a struct form.
    receivers_positions_rad = structfun(@(p) struct( ...
        'Latitude',  dms_str_2_rad(p.Latitude), ...
        'Longitude', dms_str_2_rad(p.Longitude), ...
        'Height',    p.Height), receivers_positions_dms, 'UniformOutput', false);

    training_data_geom_params.rx_pos_dms = receivers_positions_dms;
    training_data_geom_params.rx_pos_rad = receivers_positions_rad;
    % `rx_vel` is configured to be a array with all elements equal to zero,
    % given that we'll focus on ionospheric monitoring stations on GPS
    % Solutions paper.
    training_data_geom_params.rx_vel = [0; 0; 0];
    training_data_geom_params.date_time = [2014 01 02 10 00 00]; % Date/time as [YYYY MM DD hh mm ss]; consider converting to datetime in future versions.
    training_data_geom_params.simulation_time = 300;
    training_data_geom_params.dt = 0.01; % Sampling interval;
    training_data_geom_params.ipp_height = 350000;
    drift_velocities = linspace(25,125,drift_velocities_amount);
    training_data_geom_params.drift_velocity = [zeros(1,drift_velocities_amount); drift_velocities; zeros(1,drift_velocities_amount)];
end