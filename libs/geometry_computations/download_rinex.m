function rinex = download_rinex(cspsm_root_dir, datetime)
% download_rinex Returns the RINEX file as a structure.
%
% Syntax:)
%   rinex = download_rinex(dt)
%
% Description:
%   Returns a struct containing the downloaded RINEX file information.
%   Each field is a constellation, and the values is a timetable with 
%   the values for each ephemerides. This function tries to download the
%   corresponding RINEX file from https://cddis.nasa.gov/archive/gnss/data/daily.
%   The downloaded RINEX files are always stored at `cspsm_root_dir/cache`,
%   where `cspsm_root_dir` is the root directory of `cspsm.m` on you
%   system.
%
% Inputs:
%   cspsm_root_dir   - TODO:
%
%   datetime         - TODO:
%
% Notes:
%   - Downloading RINEX files from CDDIS requires user authentication. The
%     function will prompt for a valid username and password; there is no
%     anonymous or unauthenticated access.
%
% Author:
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com


%% initialization
% SEE: https://www.ordnancesurvey.co.uk/documents/resources/rinex-file-naming.pdf

cddis_url_prefix = 'https://cddis.nasa.gov/archive/gnss/data/daily/';
% day of year in the format DDD
ddd = sprintf('%03s', num2str(day(datetime, 'dayofyear')));
% data source
% NOTE: apparently, From 2020 onwards, the unique data source of the RINEX
% files was the stations. Up to 2018, the unique data source of the RINEX
% files was the receivers. 2019 was the unique year that had both receivers
% and stations as data sources.
if year(datetime) < 2020
    data_source = 'R'; % receiver
else
    data_source = 'S'; % station
end
% year in the format YYYY
yyyy = num2str(year(datetime));
% hour and minture in the format HH and MM
% NOTE: all navigation RINEX files are broadcast at 00h00m
hh = '00';
mm = '00';
% RINEX file frequency
% NOTE: RINEX navigation files seem to be broadcast daily
rinex_freq = '01D';
% contellation
constellation = 'M'; % mixed constellation
% data type indicator: N = navigation
data_type = 'N';

% file name
filename = ['BRDM00DLR_' data_source '_' yyyy ddd hh mm '_' ...
    rinex_freq '_' constellation data_type '.rnx'];

% CDDIS URL
cddis_url = [cddis_url_prefix yyyy '/brdc/' filename '.gz'];

full_filepath = fullfile(cspsm_root_dir, 'cache' , filename);

%% Download RINEX file
disp(['If you do not believe us, download it ' ...
    'manually and rerun with the RINEX file path.']);
username = input('Write the username: ', 's');
disp(['You must have wget installed in order to download the ' ...
    'RINEX file. On Linux and MacOS, it is tipically included. ' ...
    'On Windows, you will need external tools such as Chocolatey' ...
    'To have it installed.'
    ])
password = input('Write the password: ', 's');
% download
system(['wget --auth-no-challenge --user=' username ' --password=' password ...
        ' -O ' full_filepath '.gz ' cddis_url]);
% extract the .gziped file
gunzip([full_filepath '.gz']);
% delete the .gziped file
delete([full_filepath '.gz']);
% return rinex file
rinex = rinexread(filename);
end

