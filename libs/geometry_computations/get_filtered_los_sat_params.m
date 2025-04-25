function [filtered_los_sat_params, sim_params_constellations, sim_params_freqs] = get_filtered_los_sat_params(log, los_sat_params, los_prns, input_constellations, input_freq_names, trange)
%GET_LOS_SATS Summary of this function goes here
%   Detailed explanation goes here

%% Get available contellations
% all possible constellations and their repectives IDS, shown in
% `los_sat_params.Source`
all_constellations = ["gps","galileo"]; % "glonass", "beidou"
all_ids = ["PRN:","GAL Sat ID:"];

% get all available contellations shown in `los_sat_params.Source`
available_constellations = all_constellations(arrayfun(@(id) any(contains(los_sat_params.Source, id)), all_ids ));

if isempty(char(available_constellations))
    log.error('', ['For the given prograpation geometry, datetime, and RINEX file, there\n' ...
        'is no in-view satellites, at all. Try a different combination of input arguments\n' ...
        'to find visible satellites.'])
end

%% resolve simulation's contellation parameter: from either "all", user argin, or interective choise based on available constellations in the RINEX file

% if the user passed no constellation, they must decide on that at runtime
% interactive constellation selection
if isempty(char(input_constellations))
    fprintf(['You have not passed any constellation. For this RINEX file,\n'...
        'datetime, and propagation geometry, the following are availables:\n\n']);

    for i = 1:numel(available_constellations)
        fprintf('%d) %s\n', i, available_constellations(i));
    end
    select = input(['\nPass a comma-separated number sequence ' ...
        '(no spaces in between) to select which constellation(s) ' ...
        'you want\n'], "s");
    
    % get the selected indices
    select = cellfun(@(x) double(string(x)), strsplit(select, ','));
    % create a string array with the selected constellation(s)
    sim_params_constellations = available_constellations(select);
elseif any(contains(string(input_constellations), ["glonass", "beidou"]))
    log.error('', ['Sorry, but at the moment the GLONASS and Beidous are not implemented\n' ...
                    'for the satellite() function. Try either GPS or Galileo. See more in\n' ...
                    'https://www.mathworks.com/help/satcom/ref/satellitescenario.satellite.html']);
elseif input_constellations == "all"
    sim_params_constellations = available_constellations;
elseif ~any(ismember(input_constellations, available_constellations))
    log.error('', ['You have chosen constellation(s) that are not available for the RINEX file,\n' ...
        'datetime, and propagation geometry. Try to use a different conbination of input parameters\n' ...
        'to match the input with available contellations.'])
else
    assert(any(ismember(input_constellations, all_constellations)), 'The user input %s should be member of all possible contellations: %s', join(input_constellations, ','), join(all_constellations, ','));
    sim_params_constellations = input_constellations;
end

%% resolve simulation's frequency parameters: from either "all", user argin
% frequency names
all_freq_names.gps     = {'L1', 'L2', 'L5'};
all_freq_names.galileo = {'E1', 'E5a', 'E5b', 'E6'};
all_freq_names.glonass = {'G1', 'G2', 'G3'};
all_freq_names.beidou  = {'B1', 'B2', 'B3'};

% frequency values
% GPS L1 (1575.42e6), L2 (1227.60e6), L5 (1176.45e6) relative to fundamental
gps_f0 = 10.23e6; % GPS fundamental frequency
all_freq_values.gps     = [154*gps_f0, 120*gps_f0, 115*gps_f0];
% Galileo E1, E5a, E5b, E6
all_freq_values.galileo = [1575.42e6, 1176.45e6, 1207.14e6, 1278.75e6];
% GLONASS G1, G2, G3
all_freq_values.glonass = [1602e6, 1246e6, 1202.025e6];
% BeiDou B1, B2, B3
all_freq_values.beidou  = [1561.098e6,1207.14e6, 1268.52e6];

% create a container me: freq names -> freq values
all_freq_names.all_constellations = [all_freq_names.gps all_freq_names.galileo ...
    all_freq_names.glonass all_freq_names.beidou];
all_freq_values.all_constellations = [all_freq_values.gps all_freq_values.galileo ...
    all_freq_values.glonass all_freq_values.beidou];
name2value_freq = containers.Map(all_freq_names.all_constellations, all_freq_values.all_constellations);


% if input frequency name is `'all'`, replace it by all frequencies from
% the simulation paramater's constellations
if strcmpi(input_freq_names, 'all')
    input_freq_names = all_freq_names.all_constellations;
else

end

% select frequencies
for sim_params_constellation=string(sim_params_constellations)
    sim_params_freqs.(sim_params_constellation) = [];
    
    for freq_name_for_this_const=input_freq_names(ismember(input_freq_names, all_freq_names.(sim_params_constellation)))
        sim_params_freqs.(sim_params_constellation) = [sim_params_freqs.(sim_params_constellation) name2value_freq(freq_name_for_this_const{:})];
    end
end

%% Get constellation-filtered LOS satellites
filtered_ids = all_ids(contains(all_constellations, sim_params_constellations));

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

%% Get PRN-filtered LOS satellites
% TODO: also filter by PRN using `los_prns`

end

