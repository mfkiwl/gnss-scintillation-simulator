function rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_parameters)
% get_rhof_veff_ratio Computes the scaling parameter (ρ_F/v_eff) for scintillation simulation.
%
%   rhoF_veff_ratio = get_rhof_veff_ratio(general_parameters) returns the mean
%   ratio between the square root of the effective ionospheric piercing point (IPP)
%   range and the effective IPP velocity (v_eff) scaled by sqrt(2*pi*L1 frequency). 
%   This scaling parameter is used in scintillation simulation as described in [1].
%
% Notes: 
% - I've tried to maintain a similar methodology here from the code
%   `RunPropGeomCalc.m`, who main author was Charles Rino 
%   (https://www.researchgate.net/profile/Charles-Rino), and it was 
%   modified by Dongyang Xu (https://www.researchgate.net/profile/Noah-Xu) 
%   and Yu Jiao (https://www.researchgate.net/profile/Yu-Jiao-2/research).
%   There are still many improvements that could be achieved here.
%
% Dependencies:
%   - GenUserTraj
%   - ExtractRINEXeph
%   - UT2GPStime (a function written by Charles Rino; TODO (Rodrigo): Consider substituting this later, since this seems to be very old.)
%   - PropGeomCalc
%
% References:
% [1] Jiao, Yu, Rino, Charles, Morton, Yu (Jade), Carrano, Charles, 
%     "Scintillation Simulation on Equatorial GPS Signals for Dynamic 
%     Platforms," Proceedings of the 30th International Technical Meeting 
%     of the Satellite Division of The Institute of Navigation (ION GNSS+ 
%     2017), Portland, Oregon, September 2017, pp. 1644-1657. 
%     https://doi.org/10.33012/2017.15258
%
% See also: UT2GPStime, GenUserTraj, ExtractRINEXeph, PropGeomCalc
%
% Written by: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    %% Get the user trajectory (receiver's LLH coordinates)
    origin_llh = GenUserTraj(general_parameters);
    
    %% Extract the ephemeris for the GPS satellites
    eph = ExtractRINEXeph(general_parameters);

    %% Generate 1-second time samples for the propagation geometry calculation
    % UT2GPStime is a function written by Charles Rino.
    % TODO (Rodrigo): Consider substituting UT2GPStime later, since it appears to be very old.
    [gps_time_sec_start, gps_week, ~, ~] = UT2GPStime(general_parameters.date_time);
    gps_time_sec_end = gps_time_sec_start + general_parameters.simulation_time - 1;
    
    % Create an array with 2 rows and simulation_time columns:
    %   - The first row contains the GPS week numbers.
    %   - The second row contains the seconds counting from the start of the week.
    gps_week_in_seconds(1, :) = ones(1, general_parameters.simulation_time) * gps_week;
    gps_week_in_seconds(2, :) = gps_time_sec_start : gps_time_sec_end;
    
    % Handle the case when a GPS week is finishing (i.e., seconds exceed 604800)
    index_new_week = find(gps_week_in_seconds(2, :) >= 604800);
    % Correct the week number for these samples.
    gps_week_in_seconds(1, index_new_week) = gps_week_in_seconds(1, index_new_week) + 1;
    % Adjust the seconds to wrap around.
    gps_week_in_seconds(2, index_new_week) = gps_week_in_seconds(2, index_new_week) - 604800;

    %% Compute the propagation geometry and IPP parameters
    sat_geom = PropGeomCalc( ...
        gps_week_in_seconds, ...
        general_parameters.date_time, ...
        eph, ...
        origin_llh, ...
        general_parameters.ipp_height, ...
        general_parameters.rx_vel, ...
        general_parameters.drift_velocity);
    
    sat_receiver_range = sat_geom.sat_rnge;
    ipp_range = sat_geom.rngp;
    veff = sat_geom.veff;

    % The effective IPP range is obtained by the equation (13) of [1]:
    effective_ipp_range = ipp_range .* (sat_receiver_range - ipp_range) ./ sat_receiver_range;

    % The scaling parameter ρ_F / v_eff for the L1 frequency band is 
    % obtained by the equation (12) of [1]. Note that it is feasible to 
    % take the mean of all samples of ρ_F / v_eff, since it doesn't change
    % significantly over a few minutes. The usage of the mean value was 
    % first introduced by "Joy", in line 49 of "RunGenScintFieldRelization"
    % from the `gnss-scintillation-simulator-2-param` repository.
    rhof_veff_ratio = mean( ...
        sqrt(effective_ipp_range) ./ (veff * sqrt(2*pi * general_parameters.gps_bands(1))) ...
    );
end
