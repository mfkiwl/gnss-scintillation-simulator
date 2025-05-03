function rhof_veff_ratio_ref = get_scaling_param(rx, sat, sim_params)
%get_scaling_param Compute reference scaling parameter (rho_F/v_eff) at L1.
%
%   rhof_veff_ratio_ref = get_scaling_param(rx, sat, sim_params)
%   calculates the mean ratio of Fresnel scale to effective velocity at the
%   model's reference frequency using IPP geometry and drift velocity.
%
% Inputs:
%   rx                     - (required, platform object) Receiver platform.
%   sat                    - (required, satellite object) Satellite platform.
%   sim_params             - (required, struct) Simulation parameters with fields:
%                             .cte.spectral.freq_ref.value (Hz),
%                             .cte.c (speed of light, m/s),
%                             .ipp_altitude (m),
%                             .drift_vel_ned (m/s, 1x3 array).
%
% Output:
%   rhof_veff_ratio_ref    - (scalar) Mean (rho_F/v_eff) ratio at reference frequency.
%
% Reference:
%   [1] Jiao, Y., Rino, C., Morton, Y. (Jade), Carrano, C., “Scintillation Simulation
%       on Equatorial GPS Signals for Dynamic Platforms,” Proceedings of the 30th
%       International Technical Meeting of the Satellite Division of The Institute of
%       Navigation (ION GNSS+ 2017), Portland, OR, September 2017, pp. 1644–1657.
%       https://doi.org/10.33012/2017.15258
%
% Author:
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com
%
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

%% Initalization
% get receiver LLA (latitude [deg], longitude [deg], altitude [m]) trajectory and UTC time
[rx_traj_lla, rx_vel_ned, time_utc] = states(rx, 'CoordinateFrame','geographic');
% get satellite LLA (latitude [deg], longitude [deg], altitude [m]) trajectory
[sat_traj_lla, sat_vel_ned, ~]= states(sat, 'CoordinateFrame','geographic');
% drift velocity
drift_vel_ned = sim_params.drift_vel_ned;
% IPP altitude
ipp_altitude = sim_params.ipp_altitude;
% reference frequency for which we are computing ρF/veff
freq_ref = sim_params.cte.spectral.freq_ref.value;
% spped of light
c = sim_params.cte.c;

% convert LLA from deg to rad
rx_traj_lla(1:2,:)  = rx_traj_lla(1:2,:)*pi/180;
sat_traj_lla(1:2,:) = sat_traj_lla(1:2,:)*pi/180;

%% Compute the propagation geometry and IPP parameters
[rx2sat_range, ipp2rx_range, veff] = get_veff( ...
    time_utc, ...
    rx_traj_lla, ...
    rx_vel_ned, ...
    sat_traj_lla, ...
    sat_vel_ned, ...
    ipp_altitude, ...
    drift_vel_ned);


% Effective IPP range [1, Eq. (13)]
effective_ipp_range = ipp2rx_range .* (rx2sat_range - ipp2rx_range) ./ rx2sat_range;

% Mean ratio (rho_F / v_eff) at L1 [1, Eq. (12)]
% NOTE: The usage of the average ratio follows the approach introduced by
% "Joy" in the gnss-scintillation-simulator-2-param repository.
rhof_veff_ratio_ref = mean( ...
    sqrt(effective_ipp_range) ./ (veff * sqrt((2 * pi * freq_ref) / c)) ...
    );

end
