function rhof_veff_ratio_L1 = get_rhof_veff_ratio(gen_params)
% get_rhof_veff_ratio
%
% Syntax:
%   rhof_veff_ratio_L1 = get_rhof_veff_ratio(gen_params)
%
% Description:
%   Computes the mean value of the scaling parameter (rho_F / v_eff) for the L1
%   frequency band. This parameter is used in scintillation simulation to relate
%   the effective ionospheric piercing point (IPP) range and velocity to the
%   L1 carrier frequency. Under the hood, this function:
%       1) Generates user trajectory information (GenUserTraj).
%       2) Extracts ephemeris data (ExtractRINEXeph).
%       3) Converts UTC date/time to GPS time (UT2GPStime).
%       4) Computes satellite-to-receiver geometry (PropGeomCalc), including 
%          the IPP range and the effective IPP velocity (veff).
%       5) Applies equations (12) and (13) from [1] to obtain sqrt(effective_ipp_range) / (veff * sqrt(2*pi*f_L1/c)).
%          This ratio is averaged over the entire simulation interval and returned as
%          rhof_veff_ratio_L1.
%
% Inputs:
%   gen_params - Struct containing the required parameters and settings for
%                trajectory and geometry calculations. Common fields include:
%       .date_time      : MATLAB date vector [year, month, day, hour, minute, second]
%       .prn            : Satellite PRN of interest
%       .simulation_time: Duration for which data is computed (seconds)
%       .gps_bands      : [L1, L2, L5] frequencies in Hz
%       .c              : Speed of light (m/s)
%       .ipp_height     : Ionospheric pierce point height (m)
%       .rx_pos         : Receiver position in [lat [rad], lon [rad], height [m])
%       .rx_vel         : Receiver velocity (m/s)
%       .drift_velocity : Ionospheric drift velocity (m/s)
%       (additional fields may be present for GenUserTraj, ExtractRINEXeph, etc.)
%
% Outputs:
%   rhof_veff_ratio_L1 - Scalar representing the mean ratio (rho_F / v_eff) at L1,
%                        used for scaling in scintillation simulation.
%
% Notes:
%   - This function relies on external programs for orbit propagation and user
%     trajectory (GenUserTraj, ExtractRINEXeph) as well as geometry and IPP
%     calculations (PropGeomCalc).
%   - The conversion from UTC to GPS time is handled by UT2GPStime (originally
%     written by Charles Rino). This code section may be updated if a more
%     modern approach becomes available.
%   - The usage of an average value for (rho_F / v_eff) is an approach introduced 
%     in "RunGenScintFieldRelization" from the 
%     gnss-scintillation-simulator-2-param repository.
%
% Dependencies:
%   - UT2GPStime      : Converts UTC date/time to GPS week and seconds-of-week.
%   - GenUserTraj     : Generates user (receiver) trajectory data.
%   - ExtractRINEXeph : Extracts GPS ephemeris from RINEX files.
%   - PropGeomCalc    : Computes satellite geometry and IPP parameters.
%
% References:
%   [1] Jiao, Yu, Rino, Charles, Morton, Yu (Jade), Carrano, Charles, 
%       "Scintillation Simulation on Equatorial GPS Signals for Dynamic 
%       Platforms," Proceedings of the 30th International Technical Meeting 
%       of The Institute of Navigation (ION GNSS+ 2017), Portland, Oregon, 
%       September 2017, pp. 1644-1657. https://doi.org/10.33012/2017.15258
%
% See also:
%   UT2GPStime, GenUserTraj, ExtractRINEXeph, PropGeomCalc
%
% Example:
%   % Example usage:
%   gen_params.date_time       = [2023, 01, 10, 12, 00, 00];
%   gen_params.prn             = 18; % Satellite prn
%   gen_params.simulation_time = 300;             % 5 minutes
%   gen_params.gps_bands       = [1.57542e9, 1.22760e9, 1.17645e9];
%   gen_params.c               = 3e8;             % Speed of light
%   gen_params.ipp_height      = 350e3;           % 350 km
%   gen_params.rx_vel          = [0, 0, 0];       % Receiver at rest
%   gen_params.drift_velocity  = [0, 50, 0];         % 50 m/s eastward drift
%   ratio_L1 = get_rhof_veff_ratio(gen_params);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    %% Get the user trajectory (receiver's LLH coordinates)
    origin_llh = GenUserTraj(gen_params);
    
    %% Extract the ephemeris for the GPS satellites
    eph = ExtractRINEXeph(gen_params);

    %% Generate 1-second time samples for the propagation geometry calculation
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(gen_params.date_time);
    gps_time_sec_end = gps_time_sec_start + gen_params.simulation_time - 1;
    
    gps_week_in_seconds(1, :) = ones(1, gen_params.simulation_time) * gps_week;
    gps_week_in_seconds(2, :) = gps_time_sec_start : gps_time_sec_end;
    
    % Handle the case when crossing into a new GPS week
    index_new_week = find(gps_week_in_seconds(2, :) >= 604800);
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;

    %% Compute the propagation geometry and IPP parameters
    sat_geom = PropGeomCalc( ...
        gps_week_in_seconds, ...
        gen_params.date_time, ...
        eph, ...
        origin_llh, ...
        gen_params.ipp_height, ...
        gen_params.rx_vel, ...
        gen_params.drift_velocity);
    
    sat_receiver_range = sat_geom.sat_rnge;
    ipp_range = sat_geom.rngp;
    veff = sat_geom.veff;

    % Effective IPP range [1, Eq. (13)]
    effective_ipp_range = ipp_range .* (sat_receiver_range - ipp_range) ./ sat_receiver_range;

    % Mean ratio (rho_F / v_eff) at L1 [1, Eq. (12)]
    % The usage of the average ratio follows the approach introduced by "Joy" 
    % in the gnss-scintillation-simulator-2-param repository.
    rhof_veff_ratio_L1 = mean( ...
        sqrt(effective_ipp_range) ./ (veff * sqrt((2 * pi * gen_params.gps_bands(1)) / gen_params.c)) ...
    );

end
