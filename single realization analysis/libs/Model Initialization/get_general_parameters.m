function general_parameters = get_general_parameters(varargin)
% get_general_parameters Returns a struct with general simulation parameters.
%
% Syntax:
%   general_parameters = get_general_parameters()
%   general_parameters = get_general_parameters('RXPos', value, 'RXVel', value, ...
%                           'date_time', value, 'PRN', value, ...
%                           'simulation_time', value, 'Dt', value)
%
% Description:
%   Returns a struct containing various simulation parameters. If the user
%   supplies values for 'RXPos', 'RXVel', 'date_time', 'PRN', 'simulation_time',
%   or 'Dt' via name-value pairs, those values will be used; otherwise, the
%   default values are applied.
%
% Inputs (optional, as name-value pairs):
%   'RXPos'           - Receiver position as a 3x1 vector:
%                           [latitude (radians); longitude (radians); height (meters)].
%                       Default: [0.3876; 1.9942; 59.6780].
%
%   'RXVel'           - Receiver velocity as a 3x1 vector [V1; V2; V3], where:
%                           V1 = east-west velocity on the earth arc (m/s, eastward +),
%                           V2 = north-south velocity on the earth arc (m/s, northward +),
%                           V3 = up-down velocity (m/s, up +).
%                       Default: [100; 0; 0].
%
%   'date_time'       - Date and time as a vector [YYYY MM DD hh mm ss] or a datetime object.
%                       Default: [2014 01 02 10 00 00].
%
%   'PRN'             - Satellite PRN number (an integer). Must be within [0, 32].
%                       Default: 18.
%
%   'simulation_time' - Total simulation time in seconds.
%                       Default: 300.
%
%   'Dt'              - Sampling time in seconds.
%                       Default: 0.01.
%
% Outputs:
%   general_parameters - A struct containing the following fields:
%       c               - Speed of light in vacuum (m/s).
%       Dt              - Sampling time (s).
%       GPS_bands       - Frequencies for GPS bands [L1, L2, L5] in Hz.
%       date_time       - Date and time for the simulation.
%       simulation_time - Total simulation time (s).
%       RXPos           - Receiver position.
%       RXVel           - Receiver velocity.
%       PRN             - Satellite PRN.
%
% Example:
%   % Use default parameters:
%   params = get_general_parameters();
%
%   % Override RXPos, PRN, simulation time, and sampling time:
%   params = get_general_parameters('RXPos', [0.5; 2.1; 60], 'PRN', 12, ...
%                                   'simulation_time', 600, 'Dt', 0.005);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    %% Define default values
    defaultRXPos          = [0.3876; 1.9942; 59.6780];  % [latitude (rad); longitude (rad); height (m)]
    defaultRXVel          = [100; 0; 0];                 % [V1, V2, V3] where: V1 = east-west, V2 = north-south, V3 = up-down velocity
    defaultDateTime       = [2014 01 02 10 00 00];       % Date/time as [YYYY MM DD hh mm ss]; consider converting to datetime in future versions.
    defaultPRN            = 18;                          % Satellite PRN (must be within [0, 32])
    defaultSimulationTime = 300;                         % Total simulation time in seconds
    defaultDt             = 0.01;                        % Sampling time in seconds

    %% Create input parser and add parameters
    p = inputParser;
    
    % Add RXPos parameter with validation: must be a numeric 3-element vector.
    addParameter(p, 'RXPos', defaultRXPos, ...
        @(x) isnumeric(x) && isvector(x) && numel(x)==3);
    
    % Add RXVel parameter with validation: must be a numeric 3-element vector.
    addParameter(p, 'RXVel', defaultRXVel, ...
        @(x) isnumeric(x) && isvector(x) && numel(x)==3);
    
    % Add date_time parameter: accepts either a numeric vector of 6 elements or a datetime object.
    addParameter(p, 'date_time', defaultDateTime, ...
        @(x) (isnumeric(x) && numel(x)==6) || isdatetime(x));
    
    % Add PRN parameter with validation: must be a numeric scalar within the range [0, 32].
    addParameter(p, 'PRN', defaultPRN, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 32);
    
    % Add simulation_time parameter: must be a positive numeric scalar.
    addParameter(p, 'simulation_time', defaultSimulationTime, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    
    % Add Dt parameter: must be a positive numeric scalar.
    addParameter(p, 'Dt', defaultDt, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    
    % Parse the inputs provided by the user.
    parse(p, varargin{:});
    
    %% Build the general_parameters struct
    general_parameters = struct();
    
    % Constant simulation parameters:
    general_parameters.c = 299792458;                  % Speed of light in vacuum (m/s)
    
    % GPS frequencies: defined as multiples of the 10.23 MHz fundamental frequency.
    general_parameters.GPS_bands = [154*10.23e6, ...    % L1 frequency in Hz
                                    120*10.23e6, ...    % L2 frequency in Hz
                                    115*10.23e6];       % L5 frequency in Hz
    
    % Parameters that may be overridden by the user:
    general_parameters.RXPos           = p.Results.RXPos;           % [latitude (rad); longitude (rad); height (m)]
    general_parameters.RXVel           = p.Results.RXVel;           % [V1, V2, V3] with: V1 = east-west, V2 = north-south, V3 = up-down velocity.
    general_parameters.date_time       = p.Results.date_time;       % Simulation start date and time.
    general_parameters.PRN             = p.Results.PRN;             % Satellite PRN (0 to 32)
    general_parameters.simulation_time = p.Results.simulation_time; % Total simulation time (s)
    general_parameters.Dt              = p.Results.Dt;              % Sampling time (s)
end
