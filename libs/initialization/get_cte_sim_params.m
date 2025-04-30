function sim_params = get_cte_sim_params()
% get_sim_params Initialize a struct with constant simulation parameters.
%
% Syntax:
%   sim_params = get_cte_sim_params()
%
% Description:
%   This function ...
%
%
% Output
%   Constant parameters:
%       c               - (m/s, scalar) Speed of light in vacuum.
%       earth_radius    - (m, scalar) Earth radius
%   User-defined (or defaulted) parameters:
%       frequency       - (string or string array) Frequencies.
%       rx              - (struct) Receiver information:
%                           .vel: rx velocity
%                             .westeast: linear west-east velocity (m/s, eastward +, scalar)
%                             .southnorth: linear south-north velocity (m/s, northward +, scalar)
%                             .downup: linear down-up velocity (m/s, upward +, scalar)
%       sim_time - (s, scalar) Total simulation time.
%       dt              - (s, scalar) Sampling time.
%       ipp_altitude      - (m, scalar) IPP altitude.
%       drift_vel       - (m/s, struct) Ionosphere drift velocity.
%                           .x: TODO
%                           .y: TODO
%                           .z: TODO
%
% Notes:
%   - Check `parse_input_args()` to see the default values.
%   - Note that not all user inputs are simulation parameters. Rather, some
%   user input args are used just to produce a simulation parameter. That
%   is why we have `parse_input_args()` and `get_sim_params()` separately.
%
% [1] Jiao, Yu, Charles Rino, Yu (Jade) Morton, and Charles Carrano. 
%     “Scintillation Simulation on Equatorial GPS Signals for Dynamic
%     Platforms,” 1644–57, 2017. https://doi.org/10.33012/2017.15258.
% [2] Xu, Dongyang, Y.T. Jade Morton, Charles L. Rino, Charles S. Carrano,
%     and Yu Jiao. “A Two-Parameter Multifrequency GPS Signal Simulator for
%     Strong Equatorial Ionospheric Scintillation: Modeling and Parameter
%     Characterization.” NAVIGATION 67, no. 1 (2020): 181–95.
%     https://doi.org/10.1002/navi.350.
% [3] Carrano, Charles S., and Charles L. Rino. “A Theory of Scintillation
%     for Two-Component Power Law Irregularity Spectra: Overview and
%     Numerical Results.” Radio Science 51, no. 6 (2016): 789–813.
%     https://doi.org/10.1002/2015RS005903.
%% Physical parameters
sim_params.cte.c = 299792458;            % Speed of light in vacuum (m/s)
sim_params.cte.earth_radius = 6378.137e3; % Earth radius (m)
%% Contellation parameters
sim_params.cte.all_constellations = ["gps","galileo","glonass","beidou"];
% NOTE: at the moment, "glonass", "beidou" are not valid constellation as
% they are not implemented in the `satellite()` function
% SEE: https://www.mathworks.com/help/satcom/ref/satellitescenario.satellite.html
% TODO: when all global constellations are available, remove this field and
% sitck with `all_constellations`
sim_params.cte.valid_constellations = ["gps","galileo"];
% all possible constellations and their repectives IDS, shown in
% `los_sat_params.Source`
sim_params.cte.all_svid_prefix = ["PRN:","GAL Sat ID:"];

% frequency
% GPS L1 (1575.42e6), L2 (1227.60e6), L5 (1176.45e6) relative to fundamental
gps_f0 = 10.23e6; % GPS fundamental frequency
sim_params.cte.all_freqs.name.gps     = ["L1", "L2", "L5"];
sim_params.cte.all_freqs.value.gps     = [154*gps_f0, 120*gps_f0, 115*gps_f0];
% Galileo E1, E5a, E5b, E6
sim_params.cte.all_freqs.name.galileo = ["E1", "E5a", "E5b", "E6"];
sim_params.cte.all_freqs.value.galileo = [1575.42e6, 1176.45e6, 1207.14e6, 1278.75e6];
% GLONASS G1, G2, G3
sim_params.cte.all_freqs.name.glonass = ["G1", "G2", "G3"];
sim_params.cte.all_freqs.value.glonass = [1602e6, 1246e6, 1202.025e6];
% BeiDou B1, B2, B3
sim_params.cte.all_freqs.value.beidou  = [1561.098e6,1207.14e6, 1268.52e6];
sim_params.cte.all_freqs.name.beidou  = ["B1", "B2", "B3"];

%% spectral(AKA irregularity) parameters (U, p₁, p₂, μ₀) used as references for frequency extrapolation
% NOTE: These parameters are part of the extrapolation technique shown in
% [1, eq.(15), (16), and (17)], where the spectral parameters for L1
% (1575.42MHz) were estemated by using the IPE (iteractive parameters
% estimation), introduced in [3] and further explained in [2, Section 4].
% In this code, while the scaling parameter (ρF/veff) is initially computed
% for L1 in `get_rhof_veff_ratio()`, its U and μ₀ are hardcoded as constant
% values. Then, ρF/veff, U, and μ₀ are extrapolated from L1 to the desired
% frequency (there is no extrapolation step for p1 and p2).

% TODO: cite ref
sim_params.cte.spectral.strong.U_ref   = 2;
sim_params.cte.spectral.strong.mu0_ref = 0.55;
sim_params.cte.spectral.strong.p1      = 2.45;
sim_params.cte.spectral.strong.p2      = 3.7;
% TODO: cite ref
sim_params.cte.spectral.moderate.U_ref = 0.4;
sim_params.cte.spectral.moderate.mu0_ref = 0.7;
sim_params.cte.spectral.moderate.p1 = 2.7;
sim_params.cte.spectral.moderate.p2 = 3.3;
% TODO: cite ref
sim_params.cte.spectral.weak.U_ref   = 0.05;
sim_params.cte.spectral.weak.mu0_ref = 1;
sim_params.cte.spectral.weak.p1      = 3;
sim_params.cte.spectral.weak.p2      = 3;
% SEE: [2, Section 4]
sim_params.cte.spectral.freq_ref.name  = sim_params.cte.all_freqs.name.gps(1);   % L1
sim_params.cte.spectral.freq_ref.value = sim_params.cte.all_freqs.value.gps(1);  % L1

end

