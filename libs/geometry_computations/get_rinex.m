function [trange, rinex] = get_rinex(cspsm_root_dir, is_download_rinex, ...
    dtime, rinex_filename, sim_time)
%GET_RINEX Summary of this function goes here
%   Detailed explanation goes here

%% get RINEX from user inputs
if is_download_rinex
    rinex = download_rinex(cspsm_root_dir, ...
                dtime);

    % NOTE: Since the `dtime` was used to download from CDDIS, the start
    % simulation start time is totally defined in this variable
else
    rinex = rinexread(fullfile(cspsm_root_dir, ...
        'cache', ...
        rinex_filename));

    % NOTE: Since the RINEX file was not downloaded from CDDIS using
    % `dtime`, only hh:mm:ss of `dtime` carries meaning. Therefore,
    % the correct start_time should be DD/MM/YYYY from the RINEX file,
    % but hh:mm:ss from the `dtime`
    y = year(rinex.GPS.Time(1));
    m = month(rinex.GPS.Time(1));
    d = day(rinex.GPS.Time(1));
    h = hour(dtime);
    mn = minute(dtime);
    s = second(dtime);
    dtime = datetime(y, m, d, h, mn, s);
end

% define simulation time range
trange.start = dtime;
trange.end = dtime + seconds(sim_time);
end

