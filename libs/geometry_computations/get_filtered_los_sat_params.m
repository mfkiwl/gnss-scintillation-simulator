function filtered_los_sat_params = get_filtered_los_sat_params(log, los_sat_params, los_prns, filtered_constellations, trange)
%GET_LOS_SATS Summary of this function goes here
%   Detailed explanation goes here

%% Get available contellations
% all possible constellations and their repectives IDS, shown in
% `los_sat_params.Source`
all_constellations = ["gps","galileo"]; % "glonass", "beidou"
all_ids = ["PRN:","GAL Sat ID:"];

% get all available contellations shown in `los_sat_params.Source`
available_constellation = all_constellations(arrayfun(@(id) any(contains(los_sat_params.Source, id)), all_ids ));

if isempty(char(available_constellation))
    log.error('', ['For the given prograpation geometry, datetime, and RINEX file, there\n' ...
        'is no in-view satellites, at all. Try a different combination of input arguments\n' ...
        'to find visible satellites.'])
end

%% handle empty, `"all"`, or non-available contellation for `filtered_contellation`

% if the user passed no constellation, they must decide on that at runtime
% interactive constellation selection
if isempty(char(filtered_constellations))
    fprintf(['You have not passed any constellation. For this RINEX file,\n'...
        'datetime, and propagation geometry, the following are availables:\n\n']);

    for i = 1:numel(available_constellation)
        fprintf('%d) %s\n', i, available_constellation(i));
    end
    select = input(['\nPass a comma-separated number sequence ' ...
        '(no spaces in between) to select which constellation(s) ' ...
        'you want\n'], "s");
    
    % get the selected indices
    select = cellfun(@(x) double(string(x)), strsplit(select, ','));
    % create a string array with the selected constellation(s)
    filtered_constellations = all_constellations(select);
elseif any(contains(string(filtered_constellations), ["glonass", "beidou"]))
    log.error('', ['Sorry, but at the moment the contellations GLONASS and Beidous are implemented\n' ...
                    'for the satellite() function. Try either GPS or Galileo. See more in\n' ...
                    'https://www.mathworks.com/help/satcom/ref/satellitescenario.satellite.html']);
elseif filtered_constellations == "all"
    filtered_constellations = all_constellations;
end

%% Get user-filtered LOS satellites
filtered_ids = all_ids(contains(all_constellations, filtered_constellations));

% TODO: also filter by PRN using `los_prns`

filtered_los_sat_params = los_sat_params(contains(los_sat_params.Source, filtered_ids), :);
trange.start.Format = 'HH:mm:ss';
trange.end.Format = 'HH:mm:ss';

if any(filtered_los_sat_params.Duration < seconds(trange.end - trange.start))
    discarted_sats = filtered_los_sat_params(filtered_los_sat_params.Duration < seconds(trange.end - trange.start), :);
    for i = height(discarted_sats)
        log.info('The satellite %s has been discated as it lost LOS within the time range [%s, %s]', discarted_sats(i,:).Source, trange.start, trange.end);
    end
    filtered_los_sat_params = filtered_los_sat_params(filtered_los_sat_params.Duration == seconds(trange.end - trange.start), :);
end

end

