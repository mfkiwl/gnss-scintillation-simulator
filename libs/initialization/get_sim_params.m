function sim_params = get_sim_params(rx_vel_ned, drift_vel_ned, ipp_altitude)
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


%% Constant simulation parameters:
sim_params.c = 299792458;            % Speed of light in vacuum (m/s)
sim_params.earth_radius = 6378.137e3; % Earth radius (m)

%% User-defined or dafaulted parameters

%% drift velocity
sim_params.drift_vel_ned = drift_vel_ned;

%% receiver velocity
% TODO: use the tajectory object as a simulation parameter and use its
% velocity instead the user input
rx.vel_ned = rx_vel_ned;
sim_params.rx = rx;

%% IPP altitude
sim_params.ipp_altitude = ipp_altitude;


end

