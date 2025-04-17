function rhof_veff_ratio_L1 = get_rhof_veff_ratio(sim_params)
% get_rhof_veff_ratio
%
% Syntax:
%   rhof_veff_ratio_L1 = get_rhof_veff_ratio(sim_params)
%
% Description:
%   Computes the mean value of the scaling parameter (rho_F / v_eff) for the L1
%   frequency band. This parameter is used in scintillation simulation to relate
%   the effective ionospheric piercing point (IPP) range and velocity to the
%   L1 carrier frequency. Under the hood, this function:
%       1) Generates user trajectory information (set_rx_traj).
%       2) Extracts ephemeris data (ExtractRINEXeph).
%       3) Converts UTC date/time to GPS time (UT2GPStime).
%       4) Computes satellite-to-receiver geometry (PropGeomCalc), including 
%          the IPP range and the effective IPP velocity (veff).
%       5) Applies equations (12) and (13) from [1] to obtain sqrt(effective_ipp_range) / (veff * sqrt(2*pi*f_L1/c)).
%          This ratio is averaged over the entire simulation interval and returned as
%          rhof_veff_ratio_L1.
%
% Inputs:
%   sim_params - Struct containing the required parameters and settings for
%                trajectory and geometry calculations. Common fields include:
%       .datetime      : MATLAB date vector [year, month, day, hour, minute, second]
%       .prn            : Satellite PRN of interest
%       .sim_time: Duration for which data is computed (seconds)
%       .gps_bands      : [L1, L2, L5] frequencies in Hz
%       .c              : Speed of light (m/s)
%       .ipp_height     : Ionospheric pierce point height (m)
%       .rx_pos         : Receiver position in [lat [rad], long [rad], height [m])
%       .rx_vel         : Receiver velocity (m/s)
%       .drift_velocity : Ionospheric drift velocity (m/s)
%       (additional fields may be present for set_rx_traj, ExtractRINEXeph, etc.)
%
% Outputs:
%   rhof_veff_ratio_L1 - Scalar representing the mean ratio (rho_F / v_eff) at L1,
%                        used for scaling in scintillation simulation.
%
% Notes:
%   - This function relies on external programs for orbit propagation and user
%     trajectory (set_rx_traj, ExtractRINEXeph) as well as geometry and IPP
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
%   UT2GPStime, set_rx_traj, ExtractRINEXeph, PropGeomCalc
%
% Example:
%   % Example usage:
%   sim_params.datetime       = [2023, 01, 10, 12, 00, 00];
%   sim_params.prn             = 18; % Satellite prn
%   sim_params.sim_time = 300;             % 5 minutes
%   sim_params.gps_bands       = [1.57542e9, 1.22760e9, 1.17645e9];
%   sim_params.c               = 3e8;             % Speed of light
%   sim_params.ipp_height      = 350e3;           % 350 km
%   sim_params.rx_vel          = [0, 0, 0];       % Receiver at rest
%   sim_params.drift_velocity  = [0, 50, 0];         % 50 m/s eastward drift
%   ratio_L1 = get_rhof_veff_ratio(sim_params);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com
    
    %% Extract the ephemeris for the GPS satellites
    eph = ExtractRINEXeph(sim_params);

    %% Generate 1-second time samples for the propagation geometry calculation
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(sim_params.datetime);
    gps_time_sec_end = gps_time_sec_start + sim_params.sim_time - 1;
    
    gps_week_in_seconds(1, :) = ones(1, sim_params.sim_time) * gps_week;
    gps_week_in_seconds(2, :) = gps_time_sec_start : gps_time_sec_end;
    
    % Handle the case when crossing into a new GPS week
    index_new_week = find(gps_week_in_seconds(2, :) >= 604800);
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;

    %% Compute the propagation geometry and IPP parameters
    sat_geom = PropGeomCalc( ...
        gps_week_in_seconds, ...
        sim_params.datetime, ...
        eph, ...
        rx_traj_llh, ...
        sim_params.ipp_height, ...
        sim_params.rx_vel, ...
        sim_params.drift_velocity);
    
    sat_receiver_range = sat_geom.sat_rnge;
    ipp_range = sat_geom.rngp;
    veff = sat_geom.veff;

    % Effective IPP range [1, Eq. (13)]
    effective_ipp_range = ipp_range .* (sat_receiver_range - ipp_range) ./ sat_receiver_range;

    % Mean ratio (rho_F / v_eff) at L1 [1, Eq. (12)]
    % The usage of the average ratio follows the approach introduced by "Joy" 
    % in the gnss-scintillation-simulator-2-param repository.
    rhof_veff_ratio_L1 = mean( ...
        sqrt(effective_ipp_range) ./ (veff * sqrt((2 * pi * sim_params.gps_bands(1)) / sim_params.c)) ...
    );

end
