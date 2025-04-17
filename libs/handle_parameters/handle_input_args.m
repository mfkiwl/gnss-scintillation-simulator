function parsed_input_args = handle_input_args(varargin)
% handle_input_args Returns a struct with general simulation parameters.
%
% Syntax:
%   general_parameters = handle_input_args(varargin)
%
% Description:
%   Parses the raw name–value pair arguments provided by the user and
%   returns a single struct containing either the user-supplied values
%   or appropriate defaults. This function’s primary purpose is to
%   validate all inputs and ensure there are no conflicting or invalid
%   parameter combinations before the simulation begins. The input arg
%   parsing is accomplished in two steps. First, the independent arguments
%   are parsed, thus checking any inconsistency. Second, the remaining and
%   dependent arguments are analyzed; the valid values for the dependent
%   arguments depend on the independent argument values.
%
% Inputs (optional, as name-value pairs):
%   'rx_origin'       - (optional, [rad, rad, m], 3x1 array)
%                       Receiver position as [latitude; longitude; height].
%                       Default: [0.3876; 1.9942; 59.6780].
%   'rx_vel'          - (optional, m/s, 3x1 array) Receiver velocity
%                       as [v1; v2; v3], where:
%                           v1 = west-east velocity on the earth arc (eastward +),
%                           v2 = north-south velocity on the earth arc (northward +),
%                           v3 = up-down velocity (upward +).
%                       Default: [0; 0; 0].
%
%   'rinex_filename' -  (opitional, string) full file path of the RINEX
%                       file v3.04. If `'rinex_filename'` is `""`, then this
%                       program tries to download an ephemeris file from
%                       https://cddis.nasa.gov/archive/gnss/data/daily
%                       using the `'datetime'` argument to search for.
%                       Default: ""
%
%   'datetime'       -  (optional, datetime or NaN) A datetime object representing
%                       the specific date and time used as a reference
%                       point when searching for the a RINES file on CDDIS.
%                       If `'rinex_filename'` is not empty, this variable is
%                       ignored. We only support RINEX v3.04 or later,
%                       which was introduced in 2015—so the provided
%                       datetime must fall in the year 2016 or later.
%                       Default: datetime([2017 01 02 10 00 00])
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
%   'frequency'      -  (optional, string or string array) The desired
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
%   'prn'            -  (optional, string or string array) Satellite PRN.
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
%   'sim_time'        - (optional, seconds, scalar) Total simulation time.
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

%
% Example:
%   % Use default parameters:
%   params = handle_input_args();
%
%   % Override rx origin, prn, and constellation:
%   params = handle_input_args('rx_origin', [0; 5; 3], 'constellation', 'gps', 'prn', 'G12')
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
default_rx_origin       = [0.3876; 1.9942; 59.6780];        % [latitude (rad); longitude (rad); height (m)]
default_rx_vel          = [0; 0; 0];                        % [v1, v2, v3] where: v1 = west-east, v2 = north-south, v3 = up-down velocity
default_datetime        = NaN;                              % datetime
default_rinex_filename  = "";                               % RINEX file path. Empty string means that a RINEX file should be downloaded from CDDIS
default_prn             = "";                               % PRNs. Empty string means that it should be defined interactively
default_constellation   = "";                               % Constellations. Empty string means that it should be defined interactively
default_frequency       = "";                               % Default frequency. Empty string means that it should be defined interactively
default_sim_time        = 300;                              % total simulation time in seconds
default_dt              = 0.01;                             % sampling time in seconds
default_ipp_height      = 350e3;                            % IPP height in meters
default_drift_vel       = [0; 125; 0];                      % Ionosphere drift velocity [vdx, vdy, vdz] in m/s

%% Phase 1: parse the independent arg inputs
    function is_valid = validade_rinex_filename(x)
        try rinexinfo(x);
            is_valid = true;
        catch
            is_valid = false;
        end
    end

p1 = inputParser;
% p1.Unmatched will contain all key-values that are not parsed in this step
p1.KeepUnmatched = true;
% add rinex_filename parameter with validation: must be a an empty string or a
% string containing a valid RINEX file
addParameter(p1, 'rinex_filename',   default_rinex_filename, ...
    @(x) (ischar(x)||isstring(x)) && (validade_rinex_filename(x)||isempty(x)));
% Add rx_origin parameter: must be a numeric 3-element vector
addParameter(p1, 'rx_origin',   default_rx_origin, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add rx_vel parameter: must be a numeric 3-element vector
addParameter(p1, 'rx_vel',      default_rx_vel, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add sim_time parameter: must be a positive numeric scalar
addParameter(p1, 'sim_time', default_sim_time, ...
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


%% Phase 2: validate dependent input args
    function is_valid = validate_frequency(constellation, freq)
        % if no constellation, only allow empty freq
        if isempty(constellation)
            is_valid = isempty(freq);
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
        is_valid = (ischar(freq) || isstring(freq)) && any(strcmpi(freq, allowed));
    end

    function is_valid = validate_prn(constellation, prn)
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
            is_valid = false;
            return;
        end

        % when constellation is empty, the PRN can only be empty as well
        if isempty(constellation)
            is_valid = isempty(prn);
            return;
        end

        prn = strtrim(char(prn));        % convert to char and trim whitespace
        if numel(prn) < 2
            is_valid = false;
            return;
        end

        prefix = upper(prn(1));          % first letter
        numpart = prn(2:end);            % the digits
        % if the numpart (e.g., 17 in G17) is not a number, return as not valid
        if ~all(isstrprop(numpart,'digit'))
            is_valid = false;
            return;
        end
        numpart = str2double(numpart);

        % Define valid ranges
        switch lower(char(constellation))
            case 'gps'
                is_valid = (prefix=='G') && (numpart>=1 && numpart<=32);
            case 'galileo'
                is_valid = (prefix=='E') && (numpart>=1 && numpart<=36);
            case 'glonass'
                is_valid = (prefix=='R') && (numpart>=1 && numpart<=24);
            case 'beidou'
                is_valid = (prefix=='C') && (numpart>=1 && numpart<=37);
            case 'all'
                is_valid = ( (prefix=='G' && numpart<=32) || ...
                    (prefix=='E' && numpart<=36) || ...
                    (prefix=='R' && numpart<=24) || ...
                    (prefix=='C' && numpart<=37) );
            otherwise
                is_valid = false;
        end
    end

    function is_valid = validate_datetime(rinex_filename, datetime)
        if ~isdatetime(datetime)
            is_valid = false;
            return
        end

        if year(datetime) < 2016
            warning([ ...
                'You must must input a datetime with year 2016' ...
                'or later in order to download a RINEX file v3.04' ...
                ])
            is_valid = false;
            return
        end

        if ~isempty(char(rinex_filename))
            disp([ ...
                'Since you passed a RINEX file name and a datetime, ' ...
                'only the minute and hour of the datetime will be ' ...
                'used. The year, month, and day will be those ' ...
                'defined in the RINEX file.']);
        end
        is_valid = true;
    end

p2 = inputParser;
% p2.Unmatched will contain all key-values that are not parsed in this step
p2.KeepUnmatched = true;
% Add frequency parameter: valid strings depend on the set constellation
addParameter(p2, 'frequency', default_frequency, ...
    @(x) validate_frequency(p1.Results.constellation, x));
% Add prn parameter: valid strings depend on the set constellation
addParameter(p2, 'prn',       default_prn, ...
    @(x) validate_prn(p1.Results.constellation, x));
% Add datetime parameter: must be datetime object and `rinex_filename`  must be
% empty
addParameter(p2, 'datetime',    default_datetime, ...
    @(x) validate_datetime(p1.Results.rinex_filename, x));
% rebuild the remaining key-value list only from the unmatched pairs in p1
unmatched_key_value = reshape([fieldnames(p1.Unmatched) struct2cell(p1.Unmatched)]',1,[]);
parse(p2, unmatched_key_value{:});

% add p2 results
results.frequency = p2.Results.frequency;
results.prn       = p2.Results.prn;
results.datetime  = p2.Results.datetime;

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

%% Build the general_parameters struct from user input

% organize rx_origin in a structure
rx_origin.lat = results.rx_origin(1);
rx_origin.long = results.rx_origin(2);
rx_origin.height = results.rx_origin(3);

% organize drift_vel in a structure
drift_vel.x = results.drift_vel(1);
drift_vel.y = results.drift_vel(2);
drift_vel.z = results.drift_vel(3);

% organize rx_vel in a structure
rx_vel.westeast = results.rx_vel(1);
rx_vel.southnorth = results.rx_vel(2);
rx_vel.downup = results.rx_vel(3);

% Parameters that may be overridden by the user:
parsed_input_args.rx_origin       = rx_origin;                % rx origin
parsed_input_args.rx_vel          = rx_vel;                   % rx velocity
parsed_input_args.drift_vel       = drift_vel;                % ionosphere drift velocity (m/s)
parsed_input_args.prn             = results.prn;              % satellite PRNs
parsed_input_args.sim_time        = results.sim_time;         % total simulation time (s)
parsed_input_args.dt              = results.dt;               % sampling time (s)
parsed_input_args.ipp_height      = results.ipp_height;       % IPP height in meters
parsed_input_args.rinex_filename  = results.rinex_filename;   % RINEX file path
parsed_input_args.datetime        = results.datetime;         % datetime
parsed_input_args.constellation   = results.constellation;    % constellations
parsed_input_args.frequency       = results.frequency;        % frequencies

end
