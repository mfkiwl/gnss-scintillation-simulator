function ephs = extract_all_rinex_eph(date_time)
% extract_all_rinex_eph
%
% Syntax:
%   eph = extract_all_rinex_eph(date_time)
%
% Description:
%   Downloads and extracts the RINEX GPS ephemeris file for a given date,
%   then returns a struct with ephemeris records for all GPS satellites.
%   For each PRN, the record valid at the simulation start time is chosen.
%
% Inputs:
%   date_time - A vector [YYYY MM DD hh mm ss] representing the
%               simulation start time.
%
% Outputs:
%   eph - A struct with fields 'prn_xx' (e.g. 'prn_1', 'prn_2', etc.),
%         each containing a 21-element column vector of ephemeris data for
%         that satellite.
%
% Notes:
%   - The ephemeris file is downloaded from the CDDIS archive if not found.
%   - Requires external functions: rinexe, UT2GPStime, and the variable
%     general_parameters with field date_time.
%   - This function was inspired in the code `extractRINEXeph` available at 
%     -> https://github.com/cu-sense-lab/gnss-scintillation-simulator/blob/master/matlab/SatelliteOrbitComputationGPS/extractRINEXeph.m
%
% Example:
%   eph = extract_all_rinex_eph([2014 01 02 10 00 00]);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

year = date_time(1);

% Create a datetime object from the year, month, and day.
dt = datetime(date_time(1), date_time(2), date_time(3));

% Compute day of year.
day_of_year = day(dt, 'dayofyear');

% Build the RINEX filename and URL.
datadir = 'https://cddis.nasa.gov/archive/gnss/data/daily/';
YYYY = num2str(year);
DDD = num2str(day_of_year);
DDD = sprintf('%03s', DDD);  % Pad DDD to 3 characters.
YY = YYYY(3:4);
filename = ['brdc', DDD, '0.', YY, 'n.Z'];
RinexFile = [datadir, YYYY, '/', DDD, '/', YY, 'n/', filename];
[~, ephfile] = fileparts(RinexFile);

if ~exist(ephfile, 'file')
    if ~exist([ephfile, '.Z'], 'file')
        username = input('Write the username: ', 's');
        password = input('Write the password: ', 's');
        system(['wget --auth-no-challenge --user=', username, ...
            ' --password=', password, ' -O ', filename, ' ', RinexFile]);
        gunzip(RinexFile, pwd)
    end
    fprintf('Unzipping ephemeris file %s \n', ephfile)
    fprintf(['Some Matlab versions do not recognize .zip compression. ', ...
             'If error occurred, uncompress the .Z file manually.\n']);
    system(['uncompress ', ephfile, '.Z']);
    fprintf('Using ephemeris file %s \n', ephfile)
else
    fprintf('Using ephemeris file %s \n', ephfile)
end

% Extract the 21-element ephemeris records from the RINEX file.
[eph_all, ~] = rinexe(ephfile);

% Get the unique PRNs available.
unique_prns = unique(eph_all(1, :));

% Convert simulation start time to GPS time.
[GPStime_sec_start, ~, ~, ~] = UT2GPStime(date_time);

% Initialize an empty struct to hold ephemeris records for each PRN.
ephs = struct();

% Loop through each unique satellite PRN in the ephemeris data.
for i = 1:length(unique_prns)
    % Get the current PRN (satellite identifier).
    prn = unique_prns(i);
    
    % Find all column indices in eph_all where the first row (holding
    % the PRN numbers) equals the current PRN.
    idx = find(eph_all(1, :) == prn);
    
    % For the current PRN, extract the corresponding toe values.
    % The toe (time of ephemeris) is assumed to be stored in the last row.
    toe_values = eph_all(end, idx);
    
    % Find the index (within toe_values) of the last record whose toe is
    % less than or equal to the simulation start time (GPStime_sec_start).
    % This gives the most recent valid ephemeris for the current PRN.
    valid_idx = find(toe_values <= GPStime_sec_start, 1, 'last');
    
    % If no valid toe is found (i.e., none of the records occur before
    % the simulation start time), default to the first available record.
    if isempty(valid_idx)
        valid_idx = 1;
    end
    
    % Map the index within toe_values back to the corresponding column in
    % the original ephemeris matrix eph_all.
    col_idx = idx(valid_idx);
    
    % Construct a field name for the current PRN in the format 'PRNxx',
    % where xx is a two-digit representation of the PRN number.
    field_name = sprintf('prn_%d', prn);
    
    % Store the entire ephemeris record (all 21 elements) in the output
    % struct under the constructed field name.
    ephs.(field_name) = eph_all(:, col_idx);
end


return
