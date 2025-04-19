function rinex = filter_out_rinex(rinex, prn, constellation, frequency, trange)
%FILTER_OUT_RINEX Summary of this function goes here
%   Detailed explanation goes here

%% Filter out undesired constellations
% if the user passed no constellation, they must decide on that at runtime
if isempty(char(constellation))
    fprintf(['You have not passed any constellation. ' ...
        'For this RINEX file, the following are availables:\n\n']);

    rinex_fieldnames = fieldnames(rinex);
    for i = 1:numel(rinex_fieldnames)
        fprintf('%d) %s\n', i, rinex_fieldnames{i});
    end
    select = input(['\nPass a comma-separated number sequence ' ...
        '(no spaces in between) to select which constellation(s) ' ...
        'you want\n'], "s");
    
    % get the selected indices
    select = cellfun(@(x) double(string(x)), strsplit(select, ','));
    % create a string array with the selected constellation(s)
    constellation = cellfun(@(x) lower(string(x)), ...
        rinex_fieldnames(select));
end

% ensure that constellation is a string and not a char array
constellation = string(constellation);

rinex_fieldnames = cellfun(@string, fieldnames(rinex));
fields_to_remove = ~ismember(lower(rinex_fieldnames), constellation);
if all(fields_to_remove)
    error('There is no constellation%s %s in this RINEX file.', ...
        plural(numel(constellation)), strjoin(constellation, ", "));
end
% remove undesired constellation fields
rinex = rmfield(rinex, rinex_fieldnames(fields_to_remove));

%% Filter out undesired PRNs
if ~isempty(char(prn))
    % TODO: filter out PRNs
end

%% Filter out undesired frequencies


if ~isempty(char(frequency))
    % TODO: filter out frequencies
end

%% select the valid ephemeris range for this simulation
% filter out ephemerides out of the simulation range
trange = timerange(trange.start, trange.end);
rinex_fieldnames = fieldnames(rinex);
for i = numel(rinex_fieldnames)
    this_fieldname = rinex_fieldnames{i};
    rinex.(this_fieldname) = rinex.(this_fieldname)(trange,:);
    if numel(rinex.(this_fieldname)) == 0
        warning(['No %s ephemerides for the simulation range ' ...
            '[%s, %s].  This constellation was removed ' ...
            'from RINEX file.'], ...
            this_fieldname, start_time, end_time);
        rinex = rmfield(rinex, this_fieldname);
    else
        fprintf(['For the constellation %s, %d satellite%s will\n' ...
            'be used, each of them using %d eph%s\n]n'], ...
            this_fieldname, ...
            numel(unique(rinex.(this_fieldname).SatelliteID)), ...
            plural(numel(unique(rinex.(this_fieldname).SatelliteID))), ...
            numel(unique(rinex.(this_fieldname).Time)), ...
            plural(numel(unique(rinex.(this_fieldname).Time))));
    end
end

if numel(fieldnames(rinex)) == 0
    error(['There are no ephemerides for the considered ' ...
        'constellation, simulation range, and RINEX file inputs. ' ...
        'Try to change one of these in order to find a ' ...
        'combination that provides ephemerides.']);
end

end