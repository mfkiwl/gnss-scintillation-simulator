function general_parameters = get_general_parameters(varargin)
% get_general_parameters Returns a struct with general simulation parameters.
%
% Syntax:
%   general_parameters = get_general_parameters()
%   general_parameters = get_general_parameters('rx_origin', value, 'rx_vel', value, ...
%                           'datetime', value, 'prn', value, ...
%                           'simulation_time', value, 'dt', value, ...
%                           'ipp_height', value, 'drift_velocity', value)
%
% Description:
%   Returns a struct containing various simulation parameters. If the user
%   supplies values for 'rx_origin', 'rx_vel', 'datetime', 'prn', 'simulation_time',
%   'dt', 'ipp_height', or 'drift_velocity' via name-value pairs, those values will
%   be used; otherwise, the default values are applied.
%
% Inputs (optional, as name-value pairs):
%   'rx_origin'          - (optional, [rad, rad, m], 3x1 array)
%                       Receiver position as [latitude; longitude; height].
%                       Default: [0.3876; 1.9942; 59.6780].
%   'rx_vel'          - (optional, m/s, 3x1 array) Receiver velocity
%                       as [v1; v2; v3], where:
%                           v1 = east-west velocity on the earth arc (eastward +),
%                           v2 = north-south velocity on the earth arc (northward +),
%                           v3 = up-down velocity (upward +).
%                       Default: [0; 0; 0].
%
%   'datetime'       -  (optional, datetime) A datetime object representing
%                       the specific date and time used as a reference
%                       point when searching for the corresponding
%                       ephemeris timestamp.
%                       Default: datetime([2014 01 02 10 00 00])
%
%   'constellation'  -  (optional, string or string array) Desired
%                       constellations. Valid constellations are `"gps"`,
%                       `"galileo"`, `"glonass`, `"beidou"` for the global
%                       satellite navigation, and `"zqss"` and `"irnss"`
%                       for regional. If the user is unsure about which
%                       constellations are available for desired datetime,
%                       they must input `""`. In this case, an interactive
%                       prompt will pop up to guide the available
%                       constellation selection.
%                       Default: ""
%
%   'frequency'       - (optional, string or string array) The desired
%                       frequency. Must be specified if and only if
%                       `constellation` is not empty. Valid frequency
%                       labels for each constellation are:
%                         • GPS     : "L1", "L2", "L5"
%                         • Galileo : "E1", "E5a", "E5b", "E6"
%                         • GLONASS : "G1", "G2", "G3"
%                         • BeiDou  : "B1", "B2", "B3"
%                       The chosen frequency must be valid for the
%                       specified constellation. For example, "L1" is a
%                       valid label if `constealltion` is either
%                       `"gps"` or `"all"'`. If `frequency` is set to
%                       `"all"` then all frquency labels of the considered
%                       constellation is chosen. On the other hand, if
%                       `frequency` is set to `""`, then an interactive
%                       prompt will pop up to guide the available frequency
%                       section.
%
%
%   'prn'             - (optional, string or string array) Satellite PRN.
%                       If the user knows the exact satellites availabe for
%                       the desired datetime, they can input their PRNs.
%                       For instance, for 24-Jun-2021 14:00:00, the user
%                       may input `["G13", "G14"]` or just `"G13"`. The
%                       PRNs should be passed as an input if and only if
%                       the `constellation` is not empty. Also, the PRNs
%                       must match with the defined constellation types. In
%                       the previous example, the user must also input
%                       `"gps"` as one of the desired constellations or
%                       `"all"`. If all PRNs are wanted, the user must set
%                       the `prn` to `"all"`. If the user is unsure about
%                       which satellites are available for the desired
%                       datetime, they must input `""`. In this case, an
%                       interactive prompt will pop up to guide the
%                       available satellite selection.
%                       Default: ""
%
%   'simulation_time' - (optional, seconds, scalar) Total simulation time.
%                       Default: 300
%
%   'dt'              - (optional, seconds, scalar) Sampling time.
%                       Default: 10e-3
%
%   'ipp_height'      - (optional, meters, scalar) Ionospheric piercing
%                       point (IPP) height.
%                       Default: 350e3
%
%   'drift_velocity'  - (optional, m/s, 3x1 array) Ionosphere drift
%                       velocity as: [vdx, vdy, vdz] where:
%                         vdx: TODO:
%                         vdy: TODO:
%                         vdy: TODO:
%                       Default: [0; 100; 0].
%
% Outputs:
%   general_parameters - A struct containing the following fields:
%       c               - (m/s, scalar) Speed of light in vacuum.
%       frequency       - (string or string array) Frequencies.
%       rx              - (struct) Receiver information:
%                           .origin: rx starting position
%                             .lat: latitude (rad, scalar)
%                             .long: longlongitude (rad, scalar)
%                             .height: height (m, scalar)
%                           Default values for .origin:
%                             .lat: 0.3876
%                             .long: 1.9942
%                             .height: 59.6780
%                           .vel: rx velocity
%                             .westeast: linear west-east velocity (m/s, eastward +, scalar)
%                             .southnorth: linear south-north velocity (m/s, northward +, scalar)
%                             .downup: linear down-up velocity (m/s, upward +, scalar)
%                           Default values for .vel:
%                             .westeast: 0
%                             .southnorth: 0
%                             .downup: 0
%       datetime        - (datetime) Datetime for the simulation.
%       prn             - (string or string array) Satellite PRNs.
%       simulation_time - (s, scalar) Total simulation time.
%       dt              - (s, scalar) Sampling time.
%       ipp_height      - (m, scalar) IPP height.
%       drift_vel       - (m/s, struct) Ionosphere drift velocity.
%                           .x: TODO
%                           .y: TODO
%                           .z: TODO
%
% Example:
%   % Use default parameters:
%   params = get_general_parameters();
%
%   % Override rx origin, prn, and constellation:
%   params = get_general_parameters('rx_origin', [0; 5; 3], 'constellation', 'gps', 'prn', 'G12')
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com
%
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com

%% Define default values
default_rx_origin          = [0.3876; 1.9942; 59.6780];        % [latitude (rad); longitude (rad); height (m)]
default_rx_vel          = [0; 0; 0];                        % [v1, v2, v3] where: v1 = east-west, v2 = north-south, v3 = up-down velocity
default_datetime        = datetime([2014 01 02 10 00 00]);  % datetime
default_prn             = "";                               % PRNs. Empty string means that it should be defined interactively
default_constellation   = "";                               % Constellations. Empty string means that it should be defined interactively
default_frequency       = "";                               % Default frequency. Empty string means that it should be defined interactively
default_simulation_time = 300;                              % total simulation time in seconds
default_dt              = 0.01;                             % sampling time in seconds
default_ipp_height      = 350e3;                            % IPP height in meters
default_drift_vel       = [0; 125; 0];                      % Ionosphere drift velocity [vdx, vdy, vdz] in m/s

%% Phase 1: parse the independent arg inputs
p1 = inputParser;
% p1.Unmatched will contain all key-values that are parsed in this step
p1.KeepUnmatched = true;
% Add rx_origin parameter with validation: must be a numeric 3-element vector
addParameter(p1, 'rx_origin',   default_rx_origin, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add rx_vel parameter with validation: must be a numeric 3-element vector
addParameter(p1, 'rx_vel',      default_rx_vel, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add datetime parameter: must be datetime object
addParameter(p1, 'datetime',    default_datetime, ...
    @isdatetime);
% Add simulation_time parameter: must be a positive numeric scalar
addParameter(p1, 'simulation_time', default_simulation_time, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add dt parameter: must be a positive numeric scalar
addParameter(p1, 'dt',          default_dt, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add ipp_height parameter: must be a positive numeric scalar
addParameter(p1, 'ipp_height',  default_ipp_height, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add drift_vel parameter: must be a numeric 3-element vector.
addParameter(p1, 'drift_vel',   default_drift_vel, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% add constellation parameter: it must be one of the valid strings
addParameter(p1, 'constellation', default_constellation, ...
    @(x) (ischar(x)||isstring(x)) && (isempty(x) || any(strcmpi(x,{'gps','galileo','glonass','beidou','all'}))));
% parse
parse(p1, varargin{:});
% get results
results = p1.Results;


%% Phase 2: validate frequency & prn against the parsed constellation
    function tf = validate_frequency(constellation, freq)
        % if no constellation, only allow empty freq
        if isempty(constellation)
            tf = isempty(freq);
            return;
        end
        % define valid frequency labels for each GNSS constellation
        validFreq.gps     = {'L1', 'L2', 'L5'};
        validFreq.galileo = {'E1', 'E5a', 'E5b', 'E6'};
        validFreq.glonass = {'G1', 'G2', 'G3'};
        validFreq.beidou  = {'B1', 'B2', 'B3'};
        validFreq.always_valid = {''};
        % for "all", allow the union of the above
        validFreq.all     = unique([validFreq.gps, validFreq.galileo, validFreq.glonass, validFreq.beidou, validFreq.always_valid]);

        % convert the constellation value to lower case
        constellation = lower(char(constellation));

        % get the list of valid frequencies for the given constellation
        if isfield(validFreq, constellation)
            allowed = validFreq.(constellation);
        else
            allowed = {};
        end

        % Check that freq is a string or character vector and is among the allowed ones
        tf = (ischar(freq) || isstring(freq)) && any(strcmpi(freq, allowed));
    end

    function valid = validate_prn(constellation, prn)
        % validatePRN Returns true if `prn` is valid for the given constellation.
        %
        %  - If constellation is empty ("" or ''), only an empty PRN is allowed.
        %  - For 'gps':   PRN must be 'Gxx' where xx ∈ [1,32].
        %  - For 'galileo': PRN must be 'Exx' where xx ∈ [1,36].
        %  - For 'glonass': PRN must be 'Rxx' where xx ∈ [1,24].
        %  - For 'beidou':  PRN must be 'Cxx' where xx ∈ [1,37].
        %  - For 'all':     PRN may be any of the above.
        %
        % Ranges from EU GNSS glossary: GPS 1–32, Galileo 1–36,
        % GLONASS 1–24, BeiDou 1–37 :contentReference[oaicite:0]{index=0}.

        % Must be char or string
        if ~(ischar(prn) || isstring(prn))
            valid = false;
            return;
        end

        % when constellation is empty, the PRN can only be empty as well
        if isempty(constellation)
            valid = isempty(prn);
            return;
        end

        prn = strtrim(char(prn));        % convert to char and trim whitespace
        if numel(prn) < 2
            valid = false;
            return;
        end

        prefix = upper(prn(1));          % first letter
        numpart = prn(2:end);            % the digits
        % if the numpart (e.g., 17 in G17) is not a number, return as not valid
        if ~all(isstrprop(numpart,'digit'))
            valid = false;
            return;
        end
        numpart = str2double(numpart);

        % Define valid ranges
        switch lower(char(constellation))
            case 'gps'
                valid = (prefix=='G') && (numpart>=1 && numpart<=32);
            case 'galileo'
                valid = (prefix=='E') && (numpart>=1 && numpart<=36);
            case 'glonass'
                valid = (prefix=='R') && (numpart>=1 && numpart<=24);
            case 'beidou'
                valid = (prefix=='C') && (numpart>=1 && numpart<=37);
            case 'all'
                valid = ( (prefix=='G' && numpart<=32) || ...
                    (prefix=='E' && numpart<=36) || ...
                    (prefix=='R' && numpart<=24) || ...
                    (prefix=='C' && numpart<=37) );
            otherwise
                valid = false;
        end
    end

p2 = inputParser;
% p2.Unmatched will contain all key-values that are parsed in this step
p2.KeepUnmatched = true;
addParameter(p2, 'frequency', default_frequency, ...
    @(x) validate_frequency(p1.Results.constellation, x));
addParameter(p2, 'prn',       default_prn, ...
    @(x) validate_prn(p1.Results.constellation, x));
% rebuild the remaining key-value list only from the unmatched pairs in p1
unmatched_key_value = reshape([fieldnames(p1.Unmatched) struct2cell(p1.Unmatched)]',1,[]);
parse(p2, unmatched_key_value{:});

% add results
results.frequency = p2.Results.frequency;
results.prn       = p2.Results.prn;

%% Handle unmatched key-values inputs

% Helper to choose “s” or “” for pluralization
    function s = plural(n)
        if n>1, s = 's'; else, s = ''; end
    end

unknown_keys = fieldnames(p2.Unmatched);
if ~isempty(unknown_keys)
    % Join them into a comma‑separated list
    list = strjoin(unknown_keys, ', ');
    warning('Ignored unknown key%s: %s', ...
        plural(numel(unknown_keys)), list);
end

%% Build the general_parameters struct

% Constant simulation parameters:
general_parameters.c = 299792458;  % Speed of light in vacuum (m/s)

% Parameters that may be overridden by the user:
% general_parameters.constellation =
rx_origin.lat = results.rx_origin(1);
rx_origin.long = results.rx_origin(2);
rx_origin.height = results.rx_origin(3);
rx.pos = rx_origin;

rx_vel.westeast = results.rx_vel(1);
rx_vel.southnorth = results.rx_vel(2);
rx_vel.downup = results.rx_vel(3);
rx.vel = rx_vel;

drift_vel.x = results.drift_vel(1);
drift_vel.y = results.drift_vel(2);
drift_vel.z = results.drift_vel(3);

general_parameters.rx              = rx;                         % rx information
general_parameters.datetime        = results.datetime;            % simulation start date and time
general_parameters.prn             = results.prn;                 % satellite PRNs
general_parameters.simulation_time = results.simulation_time;     % total simulation time (s)
general_parameters.dt              = results.dt;                  % sampling time (s)
general_parameters.ipp_height      = results.ipp_height;          % IPP height in meters
general_parameters.drift_vel       = drift_vel;                  % ionosphere drift velocity (m/s)
end
