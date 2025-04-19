function [outputArg1,outputArg2] = cspsm(varargin)
% cspsm Compact scintillation phase screen model
%
% Syntax:
%   general_parameters = cspsm()
%   general_parameters = cspsm('rx_origin', value, 'rx_vel', value, ...
%                           'datetime', value, 'prn', value, ...
%                           'sim_time', value, 't_samp', value, ...
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
parsed_input_args = handle_input_args(cspsm_root_dir, varargin{:});

sim_params = get_sim_params(parsed_input_args.rx_vel, ...
    parsed_input_args.drift_vel, parsed_input_args.ipp_height);

%% get ephemerides

[time_range, rinex] = get_rinex(cspsm_root_dir, ...
    parsed_input_args.is_download_rinex, parsed_input_args.datetime, ...
    parsed_input_args.rinex_filename, parsed_input_args.sim_time);

% extract the appropriate ephemerides for the user frequency,
% constellation, PNR, and timedate settings from the RINEX file
rinex = filter_out_rinex(rinex, ...
    parsed_input_args.prn, parsed_input_args.constellation, ...
    parsed_input_args.frequency, time_range);

%% Geometry computations

% set rx trajectory
sim_params.rx = set_rx_traj(sim_params.rx, parsed_input_args.rx_origin, ...
    parsed_input_args.t_samp, time_range, sim_params.earth_radius);

end

