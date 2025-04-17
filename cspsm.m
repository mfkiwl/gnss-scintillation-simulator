function [outputArg1,outputArg2] = cspsm(varargin)
% cspsm Compact scintillation phase screen model
%
% Syntax:
%   general_parameters = cspsm()
%   general_parameters = cspsm('rx_origin', value, 'rx_vel', value, ...
%                           'datetime', value, 'prn', value, ...
%                           'sim_time', value, 'dt', value, ...
%                           'ipp_height', value, 'drift_velocity', value)
%
% Description:
%   This function is the entry point of the compact scintillation phase
%   screen model (CSPSM), which models the scintillaiton effects on GNSS
%   receivers.
%
%
%
%

%% Add to path
[cspsm_root_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cspsm_root_dir,'libs')));
addpath(genpath(fullfile(cspsm_root_dir,'cache')));

%% handle input args
parsed_input_args = handle_input_args(varargin{:});

sim_params = get_sim_params(parsed_input_args.sim_time, ...
    parsed_input_args.dt, parsed_input_args.rx_vel, ...
    parsed_input_args.drift_vel, parsed_input_args.ipp_height);

%% get RINEX from user inputs
% TODO: fill this catch header

% default RINEX file path
default_rinex_filename = sim_params.default_rinex_filename;
default_rinex_filepath = fullfile(cspsm_root_dir, 'cache', ...
    default_rinex_filename);

if isempty(char(parsed_input_args.rinex_filename))
    if isdatetime(parsed_input_args.datetime)
        % RINEX was not given but datetime was, try to download a RINEX
        % file for this datetime
        try
            sim_params.rinex = download_rinex(cspsm_root_dir, ...
        parsed_input_args.datetime);
        catch
            % download was not possible, fallback to the default RINEX file
            disp([ ...
                'Download was not possible, fallback ' ...
                'to the default RINEX file ' default_rinex_filename])
            sim_params.rinex = rinexread(default_rinex_filepath);
        end
    else
        % Neither RINEX nor a datetime was given, try read default RINEX
        % file
        try
            sim_params.rinex = rinexread(default_rinex_filepath);
        catch
            disp([ ...
                'It was not possble to read the default RINEX ' ...
                'file, try to download the default datetime.']);
            sim_params.rinex = download_rinex(cspsm_root_dir, ...
                sim_params.default_datetime);
        end
    end
else
    % a RINEX file was given, try to read it
    try
        full_filepath = fullfile(cspsm_root_dir, ...
            'cache', parsed_input_args.rinex_filename);
    catch
        if isdatetime(parsed_input_args.datetime)
            disp(['It was not possible to download the given' ...
                'RINEX file, try to download a RINEX file from' ...
                'the given datetime.']);
            sim_params.rinex = download_rinex(cspsm_root_dir, ...
                parsed_input_args.datetime);
        else
            disp(['It was not possible to download the given' ...
                'RINEX file, try to download a RINEX file from' ...
                'the default datetime.']);
            sim_params.rinex = download_rinex(cspsm_root_dir, ...
                sim_params.dafult_datetime);
        end
    end
    sim_params.rinex = rinexread(full_filepath);
end

% extract the appropriate ephemerides for the user frequency,
% constellation, PNR, and timedate settings
sim_params.rinex = filter_out_rinex(sim_params.rinex, ...
    parsed_input_args.prn, parsed_input_args.datetime, ...
    parsed_input_args.contellation);

%% Geometry computations

% set rx trajectory
sim_params.rx = set_rx_traj(get_sim_params.rx, ...
    get_sim_params.sim_time, get_sim_params.earth_radius);

end

