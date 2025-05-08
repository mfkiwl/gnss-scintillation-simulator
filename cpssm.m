function out = cpssm(varargin)
% cpssm Compact phase-screen-based scintillation model
%
% Description:
%   CPSSM (Compact phase-screen-based scintillation model) simulates
%   ionospheric scintillation effects on GNSS signals. It performs:
%     1. Input parsing of receiver origin, velocity, time span,
%        constellations, PRNs, simulation duration, sampling rate,
%        IPP altitude, drift velocity, and plotting options.
%     2. RINEX ephemeris loading and initial parameter setup.
%     3. Satellite scenario creation using satelliteScenario.
%     4. Receiver platform definition (static or dynamic).
%     5. LOS satellite determination and filtering by constellation and ID.
%     6. Spectral parameter extrapolation (U, μ₀) across GNSS frequencies.
%     7. Geometric computation of scintillation: time series, amplitude,
%        phase, detrended phase, normalized PSD, and effective path ratio
%        for each satellite-frequency pair.
%     8. Output assembly and optional visualization of scintillation realizations.
%
% Syntax:
%   out = cpssm()
%   out = cpssm('rx_origin', [lat lon alt], 'rx_vel', value, ...
%               'datetime', value, 'svid', value, ...
%               'sim_time', value, 't_samp', value, ...
%               'ipp_altitude', value, 'drift_velocity', value, ...
%               'is_plot', true)
%
% Inputs (optional, as name-value pairs):
%   'rx_origin'       - (optional, [deg, deg, m], 1x3 array)
%                       Receiver position as [latitude; longitude; altitude].
%                       Default: [0.3876; 1.9942; 59.6780].
%
%   'rx_vel_ned'      - (optional, m/s, 1x3 array) Receiver velocity
%                       as [v1; v2; v3] in NED (Noth-East-Down), that is:
%                           v1 = south-north velocity on the earth arc (northward +),
%                           v2 = west-east velocity on the earth arc (eastward +),
%                           v3 = up-down velocity (downward +).
%                       Default: [0; 0; 0].
%
%   'rinex_download' -  (optional, logical scalar) A logical value
%                       indicating whether the RINEX file should be
%                       downloaded from https://cddis.nasa.gov/gnss/data/daily.
%                       If it is true, the CPSSM downloads the RINEX file
%                       for the same DD/MM/YYYY as defined in the datetime
%                       argument.
%                       Default: false
%
%   'rinex_filename' -  (opitional, string) full file path of the RINEX
%                       file v3.04. You can pass `'rinex_filename'` if and
%                       only if `rinex_download` is set to false. If it is
%                       true, `'rinex_filename'` is ignored.
%                       Default: "BRDM00DLR_R_20170500000_01D_MN.rnx"
%
%   'datetime'       -  (optional, datetime or NaN) A datetime object 
%                       which defines the starting hh:mm:ss simulation
%                       time. The ending time is defined as the starting
%                       time plus the simulation duration. The difinition
%                       of the DD/MM/YYYY of the simulation depends on the
%                       RINEX file: if a RINEX file is downloaded, the
%                       CPSSM uses the DD/MM/YYYY defined in this datetime.
%                       Otherwise, if a local RINEX  file is used, the
%                       DD/MM/YYYY is the same as defined in the RINEX
%                       file, and the DD/MM/YYYY of this `datetime` is
%                       therefore ignored. You should note that, if you are
%                       downloading a RINEX file, you must pass a datetime
%                       whose year is 2016 or later because only RINEX
%                       v3.04 is supported.
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
%                       Default: "all"
%
%   'svid'            -  (optional, string or string array) Satellite SVID.
%                       If the user knows the exact satellites availabe for
%                       the desired datetime, they can input their PRNs.
%                       For instance, for 24-Jun-2021 14:00:00, the user
%                       may input `["G13", "G14"]` or just `"G13"`. The
%                       PRNs should be passed as an input if and only if
%                       the `constellation` is not empty. Also, the PRNs
%                       must match with the defined constellation types. In
%                       the previous example, the user must also input
%                       `"gps"` as one of the desired constellations or
%                       `"all"`. If all SVIDs are wanted, the user must set
%                       the `svid` to `"all"`. If the user is unsure about
%                       which satellites are available for the desired
%                       datetime, they must input `""`. In this case, an
%                       interactive prompt will pop up to guide the
%                       available satellite selection.
%                       Default: ""
%
%   'severity'        - (optional, scalar string): The scintillation 
%                       severity. It must be either `"weak"`, `"moderate"`,
%                       or `"strong"`. The severity affects the disturbance
%                       level caused by the scintillation realization,
%                       which in turns alters the scitilllation-related
%                       outputs.
%                       Default: "strong"
%
%   'sim_time'        - (optional, seconds, scalar) Total simulation time.
%                       Default: 300
%
%   't_samp'          - (optional, seconds, scalar) Sampling time of the
%                       scintillation time series realization.
%                       Default: 10e-3
%
%   'ipp_altitude'      - (optional, meters, scalar) Ionospheric piercing
%                       point (IPP) altitude.
%                       Default: 350e3
%
%   'drift_vel_ned'    - (optional, m/s, 1x3 array) Ionosphere drift
%                        velocity as [vdx, vdy, vdz] in NED
%                        (Noth-East-Down), that is:
%                         vdx: west-east velocity on the earth arc (eastward +),
%                         vdy: south-north velocity on the earth arc (northward +),
%                         vdy: up-down velocity (downward +).
%                       Default: [0; 100; 0].
%
%   'plot'             - (optional, logical scalar) A logical scalar
%                      indicating whether plots concerning the ionospheric
%                      scintillation realization should be shown. It does
%                      not affect the simulation.
%                      Default: false
%
%   'play'             - (optional, logical scalar) A logical scalar
%                      indicating whether a animation of the geometry
%                      between the receiver and the satellites should be
%                      shown. It does not not affect the simulation.
%                      Default: false
%
% Author:
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com
%
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

%% Initialize
% Add to path
[cpssm_root_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cpssm_root_dir,'libs')));
addpath(genpath(fullfile(cpssm_root_dir,'cache')));

[parsed_argins,sim_params,log] = initialize_cpssm(cpssm_root_dir, varargin);

%% Scenario
sim_params = set_scenario(log, cpssm_root_dir, sim_params, parsed_argins);

%% Scintillation realization
out = get_scintillation(sim_params);

%% Plot output
if parsed_argins.is_plot
    plot_scintillation_psd(cpssm_root_dir, out);
    % plot_scintillation_time_series(out);
end

%% Play satellite-receiver scenarios
if parsed_argins.is_play
    play(out.satelliteScenario);
end
end