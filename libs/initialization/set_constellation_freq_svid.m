function sim_params = set_constellation_freq_svid(log, sim_params, ...
    parsed_argins, los_sats_params)
%SET_CONSTELLATION_AND_FREQ Summary of this function goes here
%   Detailed explanation goes here
%% Initialization
all_constellations = sim_params.const.all_constellations;
all_svids = sim_params.const.all_svid_prefix;
all_freq = sim_params.const.all_freqs;
is_multiconst_sats = sim_params.is_multiconst_sats;

input_constellations = parsed_argins.constellations;
input_freq_names = parsed_argins.frequencies;
input_svids = parsed_argins.svids;

%% Get available contellations
% get all available contellations shown in `los_sat_params.Source`
available_constellations = all_constellations(arrayfun(@(id) any(contains(los_sats_params.Source, id)), all_svids ));

if isempty(char(available_constellations))
    log.error('', ['For the given prograpation geometry, datetime, and RINEX file, there\n' ...
        'is no in-view satellites, at all. Try a different combination of input arguments\n' ...
        'to find visible satellites.'])
end

%% resolve unavailability of GLONASS and Beidou

idx = ismember(input_constellations, ["glonass", "beidou", "zqss", "navic"]);
msg = ['At the moment the GLONASS, Beidou, QZSS, and NavIC are not implemented for\n' ...
    'the satellite() function. Stick with GPS and Galileo. See more in\n' ...
    'https://www.mathworks.com/help/satcom/ref/satellitescenario.satellite.html.'];

if all(idx)
    log.error('', [msg '\n\nPlease, select at least GPS or Galileo.']);
elseif any(idx)
    if is_multiconst_sats
        log.warning('', [msg '\n\nSince multiconstellation satellite mode is enabled, those GPS/Galileo satellites\n' ...
        'will be used to transmit GLONASS, Beidou, QZSS, and/or NavIC signals.']);
    else
        % remove GLONASS and Beidou
        input_constellations = input_constellations(~idx);
        log.warning('', [msg '\nFor this simulation, GLONASS, Beidou, QZSS, and/or NavIC will be disregarded.']);
        assert(numel(input_constellations) >= 1, "Unexpected empty input constellation after disregarding GLONASS, Beidou, QZSS, and/or NavIC. This is probably a logical error in the code.");
    end
end

%% resolve simulation's contellation parameter: from either "all", user argin, or interective choise based on available constellations in the RINEX file

% use all available constellations
if input_constellations == "all"
    log.info('The constellation option was set to all.\nThe following constellation are available from the RINEX file and will be used: %s', ...
        join(available_constellations, ','));
    sim_params_constellations = available_constellations;
% all input constellations are available or it is multiconstellation sat mode
elseif is_multiconst_sats || all(ismember(input_constellations, available_constellations))
    sim_params_constellations = input_constellations;
% not all input constellations are available
elseif any(ismember(input_constellations, available_constellations))
    idx = ismember(input_constellations, available_constellations);
    log.warning('', ['You have chosen the constellation(s) %s but they are not available for the RINEX file,\n' ...
        'datetime, and propagation geometry. Therefore, these constellation(s) will be dirergarded.'], ...
        join(input_constellations(~idx), ','));
    sim_params_constellations = available_constellations(idx);
% interactive constellation selection
elseif strcmpi(input_constellations, "")
    fprintf(['You have not passed any constellation. For this RINEX file,\n'...
        'datetime, and propagation geometry, the following are availables:\n\n']);

    for i = 1:numel(available_constellations)
        fprintf('%d) %s\n', i, available_constellations(i));
    end
    idx = input(['\nPass a comma-separated number sequence ' ...
        '(no spaces in between) to select which constellation(s) ' ...
        'you want\n'], "s");
    
    % get the selected indices
    idx = cellfun(@(x) double(string(x)), strsplit(idx, ','));
    % create a string array with the selected constellation(s)
    sim_params_constellations = available_constellations(idx);
    if parsed_argins.is_multiconst_sats && isscalar(sim_params_constellations)
        log.warning('', ['You have selected multiconstellation satellite mode, but have chosen only\n' ...
            'one constellation. Therefore the multiconstellation satellite mode will take no effect.'])
    end
% none of the input constellations were available
else
    assert(isempty(intersect(input_constellations, available_constellations)), ...
        'Unexpected non-empty intersection between available and input contellations at this code path.');
    
    log.error('', ['None of the requested constellation(s) [%s] are available. The available contellation(s) are: %s.\n', ...
        'Use a subset of it or try a different combination of RINEX file, datetime and propagation geometry\n' ...
        'to find another set of avaiblable contellations.'], ...
        join(input_constellations, ', '), ...
        join(available_constellations, ', '));
end

%% resolve simulation's frequency parameters: from either "all", user argin

% if input frequency name is `'all'`, replace it by all frequencies from
% the simulation paramater's constellations
if strcmpi(input_freq_names, 'all')
    input_freq_names = [all_freq.name.gps all_freq.name.galileo ...
    all_freq.name.glonass all_freq.name.beidou];
end

% select frequencies
for sim_params_constellation=sim_params_constellations
    % frequency name indices for this constellation that matches with user argin
    idx = ismember(all_freq.name.(sim_params_constellation), input_freq_names);
    if nnz(idx)
        sim_params_freq.(sim_params_constellation).name = all_freq.name.(sim_params_constellation)(idx);
        sim_params_freq.(sim_params_constellation).value = all_freq.value.(sim_params_constellation)(idx);
    else
        log.warning('', ['You have included %s as a constellation, but you have not included any %s frequency\n' ...
            'for that. Therefore, this constellation will be excluded. If you did not mean that,\n' ...
            'pass valid frequencies for that constellation.'], sim_params_constellation, sim_params_constellation);
        sim_params_constellations = sim_params_constellations(sim_params_constellations ~= sim_params_constellation);
    end
end


%% TODO: resolve SVIDS

sim_params_svids = input_svids;

%% Set values
sim_params.constellations = sim_params_constellations;
sim_params.freqs = sim_params_freq;
sim_params.svids = sim_params_svids;

end

