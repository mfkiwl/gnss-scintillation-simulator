function [sim_params_argin, download_rinex_argin] = split_input_args(varargin)
%split_input_args Split the input arguments so that it can be processed by
%the downstream function
%
% Syntax:
%   out = split_input_args(varargin)
%
% Description:
%   This funciton splits the input arguments so that it can be processed by
%   the downstream functions.
%
% Inputs (optional, as name-value pairs):
%
%   'rinex_filename' - (opitional, string) full file path of the RINEX
%                      file v3.04. If `'rinex_filename'` is `""`, then this
%                      program tries to download an ephemeris file from
%                      https://cddis.nasa.gov/archive/gnss/data/daily
%                      using the `'datetime'` argument to search for.
%                       Default: ""
%
%
%
%

% variables for the function `handle_input_args()`
sim_params_argin = {};
% variables for the function `download_rinex()`
download_rinex_argin = {};
% variables for the function `set_rx_traj()`
rx_traj_argin = {};


for i = 1:numel(varargin)/2
    switch varargin{2*i-1}
        case 'rinex_filename'
            download_rinex_argin(end+1) = varargin(2*i-1);
            download_rinex_argin(end+1) = varargin(2*i);
        case 'rx_origin'
            rx_traj_argin(end+1) = varargin(2*i-1);
            rx_traj_argin(end+1) = varargin(2*i);
        case {'rx_vel', 'prn', }
    end
end

end

