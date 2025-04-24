function [parsed_input_args, log] = parse_input_args(cspsm_root_dir, varargin)
% parse_input_args Returns a struct with general simulation parameters.
%
% Syntax:
%   general_parameters = parse_input_args(varargin)
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
%   'rx_origin'       - (optional, [deg, deg, m], 1x3 array)
%                       Receiver position as [latitude; longitude; altitude].
%                       Default: [0.3876; 1.9942; 59.6780].
%
%   'rx_vel'          - (optional, m/s, 1x3 array) Receiver velocity
%                       as [v1; v2; v3], where:
%                           v1 = west-east velocity on the earth arc (eastward +),
%                           v2 = north-south velocity on the earth arc (northward +),
%                           v3 = up-down velocity (upward +).
%                       Default: [0; 0; 0].
%
%
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
%   't_samp'              - (optional, seconds, scalar) Sampling time.
%                       Default: 1
%
%   'ipp_alt'      - (optional, meters, scalar) Ionospheric piercing
%                       point (IPP) altitude.
%                       Default: 350e3
%
%   'drift_velocity'  - (optional, m/s, 1x3 array) Ionosphere drift
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
%   params = parse_input_args();
%
%   % Override rx origin, prn, and constellation:
%   params = parse_input_args('rx_origin', [0; 5; 3], 'constellation', 'gps', 'prn', 'G12')
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
default_rx_origin       = [-23.2198 -45.8916  59.6780];         % [latitude (deg); longitude (deg); altitude (m)] -> São José dos Campos
default_rx_vel          = [0 0 0];                              % [v1, v2, v3] where: v1 = west-east, v2 = north-south, v3 = up-down velocity
default_datetime        = datetime([2021 06 24 14 00 00]);      % datetime
default_rinex_filename  = "BRDM00DLR_R_20170500000_01D_MN.rnx"; % RINEX file name.
default_is_down_rinex   = false;                                % by default, do not download a RINEX file and use either the user-defined or default RINEX file
default_prn             = "";                                   % PRNs. Empty string means that it should be defined interactively
default_constellation   = "";                                   % Constellations. Empty string means that it should be defined interactively
default_frequency       = "";                                   % Default frequency. Empty string means that it should be defined interactively
default_log_lvl         = "DEBUG";                              % Default log level
default_sim_time        = 300;                                  % total simulation time in seconds
default_t_samp          = 1;                                    % sampling time in seconds
default_ipp_alt         = 350e3;                                % IPP altitude in meters
default_drift_vel       = [0 125 0];                            % Ionosphere drift velocity [vdx, vdy, vdz] in m/s

%% Parsing phase 0: resolve the logging before anything else

p0 = inputParser;
% p.Unmatched will contain all key-values that are not parsed in this step
p0.KeepUnmatched = true;

% Add log_lvl parameter: must be either a either a nonfractional scalar, or
% a scalar string or a char array that belongs to {"DEBUG", "INFO", "WARN",
% "ERROR"}
addParameter(p0, 'log_lvl',   default_log_lvl, ...
    @(lvl) ((ischar(lvl) || (isstring(lvl) && isscalar(lvl))) && ...
    ismember(upper(string(lvl)),["DEBUG", "INFO", "WARN","ERROR"]) || ...
    (isscalar(lvl) && isnumeric(lvl) && (mod(lvl,1) == 0))) );

% parse it
parse(p0, varargin{:});

% get log
log = Logger(p0.Results.log_lvl);

% rebuilt `varargin` without `log_lvl` and its value (if they were passed)
unmatched_keys = fieldnames(p0.Unmatched);                 % e.g. {'rx_origin';'constellation';…}
unmatched_vals  = struct2cell(p0.Unmatched);                          % corresponding values
varargin = reshape([unmatched_keys unmatched_vals]',1,[]); % interleave into 1×(2*M) cell

%% Parse phase 1: parse the other the input arguments
p = inputParser;
% p.Unmatched will contain all key-values that are not parsed in this step
p.KeepUnmatched = true;
% Add download_rinex parameter: must be a logical scalar
addParameter(p, 'download_rinex',   default_is_down_rinex, ...
    @(x) isscalar(x) && islogical(x));
% Add rx_origin parameter: must be a numeric 3-element vector
addParameter(p, 'rx_origin',   default_rx_origin, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add rx_vel parameter: must be a numeric 3-element vector
addParameter(p, 'rx_vel',      default_rx_vel, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% Add sim_time parameter: must be a positive numeric scalar
addParameter(p, 'sim_time', default_sim_time, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add t_samp parameter: must be a positive numeric scalar
addParameter(p, 't_samp',          default_t_samp, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add ipp_alt parameter: must be a positive numeric scalar
addParameter(p, 'ipp_alt',  default_ipp_alt, ...
    @(x) isnumeric(x) && isscalar(x) && x>0);
% Add drift_vel parameter: must be a numeric 3-element vector.
addParameter(p, 'drift_vel',   default_drift_vel, ...
    @(x) isnumeric(x) && isvector(x) && numel(x)==3);
% add constellation parameter: it must be one of the valid strings
addParameter(p, 'constellation', default_constellation, ...
    @(x) (ischar(x)||isstring(x)) && (isempty(x) || any(strcmpi(x,{'gps','galileo','glonass','beidou','all'}))));
% add rinex_filename parameter: must be a an empty string or a
% string containing a valid RINEX file
addParameter(p, 'rinex_filename', default_rinex_filename, ...
    @(rinex_filename) validade_rinex_filename(log, cspsm_root_dir, ...
    p.Results.download_rinex, rinex_filename));
% Add frequency parameter: valid strings depend on the set constellation
addParameter(p, 'frequency', default_frequency, ...
    @(x) validate_frequency(log, p.Results.constellation, x));
% Add prn parameter: valid strings depend on the set constellation
addParameter(p, 'prn',       default_prn, ...
    @(x) validate_prn(log, p.Results.constellation, x));
% Add datetime parameter: must be datetime object and `rinex_filename`  must be
% empty
addParameter(p, 'datetime',    default_datetime, ...
    @(x) validate_datetime(log, p.Results.download_rinex, ...
    p.Results.rinex_filename, x));

% parse it
parse(p, varargin{:});

%% Handle unmatched key-values inputs
unknown_keys = fieldnames(p.Unmatched);
if ~isempty(unknown_keys)
    % Join them into a comma‑separated list
    list = strjoin(unknown_keys, ', ');
    log.warning('', 'Ignored unknown key%s: %s', ...
        plural(numel(unknown_keys)), list);
end

%% Build the general_parameters struct from user input

% Parameters that may be overridden by the user:
parsed_input_args.rx_origin         = p.Results.rx_origin;        % rx origin
parsed_input_args.rx_vel            = p.Results.rx_vel;           % rx velocity
parsed_input_args.drift_vel         = p.Results.drift_vel;        % ionosphere drift velocity (m/s)
parsed_input_args.is_download_rinex = p.Results.download_rinex;   % whether one should download the RINEX file
parsed_input_args.prn               = p.Results.prn;              % satellite PRNs
parsed_input_args.sim_time          = p.Results.sim_time;         % total simulation time (s)
parsed_input_args.t_samp            = p.Results.t_samp;           % sampling time (s)
parsed_input_args.ipp_alt           = p.Results.ipp_alt;          % IPP altitude in meters
parsed_input_args.rinex_filename    = p.Results.rinex_filename;   % RINEX file path
parsed_input_args.datetime          = p.Results.datetime;         % datetime
parsed_input_args.constellation     = p.Results.constellation;    % constellations
parsed_input_args.frequency         = p.Results.frequency;        % frequencies

end

% -------------------------------------------------------------------------

function validade_rinex_filename(log, cspsm_root_dir, ...
    is_down_rinex, rinex_filename)
if ~ischar(rinex_filename) && ~isstring(rinex_filename)
    log.error('rinex_filename must be a char or string');
end

if is_down_rinex
    log.error('', [ ...
        'You must either set `download_rinex` to true to ' ...
        'download a RINEX file from CDDIS, or provide a ' ...
        'local RINEX file name.'
        ]);
end

try
    rinexinfo(fullfile(cspsm_root_dir, 'cache', rinex_filename));
catch ME
    switch ME.identifier
        case 'nav_positioning:rinexInternal:FileNotFound'
            % Build your custom text
            log.error( ...
              ME.identifier, ...
              ['Failed to read RINEX file "%s".\n', ...
               'Please check that it exists and is readable.\n', ...
               'Original error was:\n%s'], ...
              fullfile(cspsm_root_dir, 'cache', rinex_filename), ...
              ME.message ...
            );
    end
end

end

function validate_datetime(log, is_down_rinex, rinex_filename, datetime)
if ~isdatetime(datetime)
    log.error('datetime must be a datetime object.');
end

if year(datetime) < 2016
    log.error('', [ ...
        'You must must input a datetime with year 2016' ...
        'or later in order to download a RINEX file v3.04' ...
        ])
end

if ~is_down_rinex
    log.info([ ...
        'You have passed a datetime but have not set to\n' ...
        'download a RINEX file from CDDIS. In this case, only the minute and hour\n' ...
        'of the datetime will be used. The year, month, and day will be those\n' ...
        'defined in %s. If you intended to use\n' ...
        'the datetime to download from that year, month, and day, set download_rinex\n' ...
        'to true.\n'], rinex_filename);
end
end

function validate_prn(log, constellation, prn)
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
    log.error('', 'PRN must be a string or char object.');
end

% when constellation is empty, the PRN can only be empty as well
if isempty(char(constellation)) && ~isempty(char(prn))
    log.error('', 'If you have not set a constellation, you cannot set a PRN.');
end

prn = strtrim(char(prn));        % convert to char and trim whitespace
if numel(prn) < 2
    log.error('A valid PRN must contain at least three characters.')
end

prefix = upper(prn(1));          % first letter
numpart = prn(2:end);            % the digits
% if the numpart (e.g., 17 in G17) is not a number, return as not valid
if ~all(isstrprop(numpart,'digit'))
    log.error("The two last PRN's characters must be numbers.")
end
numpart = str2double(numpart);

% Define valid ranges
switch lower(char(constellation))
    case 'gps'
        if ~((prefix=='G') && (numpart>=1 && numpart<=32))
            log.error('This PRN is not valid for GSP.')
        end
    case 'galileo'
        if ~((prefix=='E') && (numpart>=1 && numpart<=36))
            log.error('This PRN is not valid for Galileo.')
        end
    case 'glonass'
        if ~((prefix=='R') && (numpart>=1 && numpart<=24))
            log.error('This PRN is not valid for GLONASS.')
        end
    case 'beidou'
        if ~((prefix=='C') && (numpart>=1 && numpart<=37))
            log.error('This PRN is not valid for Beidou')
        end
    case 'all'
        if ~( (prefix=='G' && numpart<=32) || ...
              (prefix=='E' && numpart<=36) || ...
              (prefix=='R' && numpart<=24) || ...
              (prefix=='C' && numpart<=37) )
            log.error('There is no constellation in which this PRN is valid');
        end
end
end

function validate_frequency(log, constellation, freq)
% if no constellation, only allow empty freq
if ~ischar(freq) && ~isstring(freq)
    log.error(['You must pass a char, string, or string array ' ...
        'for the desired frequencies.'])
elseif isempty(constellation) && ~isempty(char(freq))
    log.error(['You cannot pass frequencies if you have not ' ...
        'defined a contellation'])
end
% define valid frequency labels for each GNSS constellation
validFreq.gps     = {'L1', 'L2', 'L5'};
validFreq.galileo = {'E1', 'E5a', 'E5b', 'E6'};
validFreq.glonass = {'G1', 'G2', 'G3'};
validFreq.beidou  = {'B1', 'B2', 'B3'};
validFreq.always_valid = {''};
% for "all", allow the union of the above
validFreq.all = unique([validFreq.gps, validFreq.galileo, ...
    validFreq.glonass, validFreq.beidou, validFreq.always_valid]);

% convert the constellation value to lower case
constellation = lower(char(constellation));

% get the list of valid frequencies for the given constellation
if isfield(validFreq, constellation)
    allowed = validFreq.(constellation);
else
    allowed = {};
end

% Check that freq is a string or character vector and is among the allowed ones
if ~(any(strcmpi(freq, allowed)))
    log.error(['You have passed frequency(ies) that are not valid ' ...
        'for the considered constellation(s).'])
end
end