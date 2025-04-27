function sim_params = get_cte_sim_params()
% get_sim_params Initialize a struct with general simulation parameters.
%
% Syntax:
%   sim_params = get_sim_params(rx_origin,sim_time, dt, drift_vel)
%
% Description:
%   This function initializes the struct `sim_params` from some
%   user-defined (or defaulted) input values. This structure contains all
%   simulation parameters and constants used throughout the CSPSM. While
%   this function defines simulation parameters from some user inputs (or
%   the default values), other are set during the program. Once a new field
%   is  defined, the program should not modify its value. The idea is to
%   define (only once) all parameters used be the downstream functions, and
%   store them in the `sim_params` struct.
%
% Input
%   rx_origin         - ([rad, rad, m], 3x1 array)
%                       Receiver position as [latitude; longitude; altitude].
%
%   sim_time          - (seconds, scalar) Total simulation time.
%
%   dt                - (seconds, scalar) Sampling time.
%
%   ipp_altitude        - (meters, scalar) Ionospheric piercing point (IPP)
%                       altitude.
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

%% Physical parameters
cte.c = 299792458;            % Speed of light in vacuum (m/s)
cte.earth_radius = 6378.137e3; % Earth radius (m)
%% Contellation parameters
cte.all_constellations = ["gps","galileo","glonass","beidou"];
% NOTE: at the moment, "glonass", "beidou" are not valid constellation as
% they are not implemented in the `satellite()` function
% SEE: https://www.mathworks.com/help/satcom/ref/satellitescenario.satellite.html
% TODO: when all global constellations are available, remove this field and
% sitck with `all_constellations`
cte.valid_constellations = ["gps","galileo"];
% all possible constellations and their repectives IDS, shown in
% `los_sat_params.Source`
cte.all_ids = ["PRN:","GAL Sat ID:"];

% frequency
% GPS L1 (1575.42e6), L2 (1227.60e6), L5 (1176.45e6) relative to fundamental
name.gps     = ["L1", "L2", "L5"];
gps_f0 = 10.23e6; % GPS fundamental frequency
value.gps     = [154*gps_f0, 120*gps_f0, 115*gps_f0];
% Galileo E1, E5a, E5b, E6
name.galileo = ["E1", "E5a", "E5b", "E6"];
value.galileo = [1575.42e6, 1176.45e6, 1207.14e6, 1278.75e6];
% GLONASS G1, G2, G3
name.glonass = ["G1", "G2", "G3"];
value.glonass = [1602e6, 1246e6, 1202.025e6];
% BeiDou B1, B2, B3
value.beidou  = [1561.098e6,1207.14e6, 1268.52e6];
name.beidou  = ["B1", "B2", "B3"];

all_freqs.name  = name;
all_freqs.value = value;
cte.all_freqs = all_freqs;

%% output

sim_params.cte = cte;

end

