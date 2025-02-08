function       eph=ExtractRINEXeph(general_parameters)
% USAGE:  Extract RINEX ephemeris file from website 
% https://cddis.nasa.gov/archive/gnss/data/daily/YYYY/DDD/YYn/brdcDDD0.YYn.Z
% and store the parameters in the eph struct:
%
% Modified by Rodrigo (06/02/2025):
% Substituted the `userInput` from the `gnss-scintillation-simulator-2-param`
% repository for the `general_parameters` in this modification. 
% In addition, i've substituted the usage of datenum for a more flexible 
% datetime object here.
%
% TODO: For now, the datetime object is being instantiated in this
% function. However, we should work entirely with a date time object, for
% the sake of readability and flexibility of the code. In a later version
% of the code, we should introduce this modification.

year = general_parameters.date_time(1);

% Create a datetime object from the year, month, and day components
dt = datetime(general_parameters.date_time(1), ...
              general_parameters.date_time(2), ...
              general_parameters.date_time(3));
          
% Compute the day of year directly
day_of_year = day(dt, 'dayofyear');

PRN = general_parameters.prn;

datadir = 'https://cddis.nasa.gov/archive/gnss/data/daily/';
YYYY = num2str(year);
DDD = num2str(day_of_year);
DDD = sprintf('%03s', DDD); % the format specifier %03s pads DDD with leading zeros to ensure it is at least 3 characters long

YY = YYYY(3:4);
filename = ['brdc',DDD,'0.',YY,'n.Z'];
RinexFile = [datadir,YYYY,'/',DDD,'/',YY,'n/',filename];
[~,ephfile] = fileparts(RinexFile);
if ~exist(ephfile,'file')
    if ~exist([ephfile,'.Z'],'file')
        username = input('Write the username: ', 's');
        password = input('Write the password: ', 's');

        system(['wget --auth-no-challenge --user=', username ' --password=', password, ' -O ', filename,' ', RinexFile])
        
        gunzip(RinexFile,pwd)
    end
    fprintf('Unzipping ephemeris file %s \n',ephfile)
    fprintf(['Some Matlab versions do not recognize .zip compression.' ...
    ' If error occurred, please go to the folder and manually uncompress the .Z ephemeris file.']);
    system(['uncompress ',ephfile,'.Z']);
    fprintf('Using ephemeris file %s \n',ephfile)
    %Matlab R2015a does not recognize .zip compression => manual uncompress nessary!!
else
    fprintf('Using ephemeris file %s \n',ephfile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eph,~] =rinexe(ephfile); %Extract 21XN ephemeris files
%eph(21,:) is toe for each 21-element ephemeris
neph=find(eph(1,:)==PRN);
if isempty(neph)
    error('RINEX error ')
else
    eph=eph(:,neph);
end

% truncate the eph data to right after the user input time.
[GPStime_sec_start,GPSweek,GPSweek_z,leapsec]=UT2GPStime(general_parameters.date_time);

ind = find(eph(end,:)>GPStime_sec_start,1,'first')-1;
if isempty(ind)
    ind = length(eph(end,:));
elseif ind==0
    ind = 1;
end
eph = eph(:,ind);

return