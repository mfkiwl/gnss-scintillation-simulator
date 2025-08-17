function filtered_los_sat_params = get_filtered_los_sat_params(log, sim_params, los_sats_params)
%GET_LOS_SATS Summary of this function goes here
%   Detailed explanation goes here
%% Initialization
all_constellations = sim_params.const.all_constellations;
% sim_params_svids = sim_params.svids;
% TODO: remove this when GLONASS, Beidou, ZQSS, or NavIC become available
idx = ismember(sim_params.constellations, ["glonass", "beidou", "zqss", "navic"]);
sim_params_constellations = sim_params.constellations(~idx);
all_svid_prefix = sim_params.const.all_svid_prefix;
time_start = sim_params.temporal_support(1);
time_end = sim_params.temporal_support(end);

%% Get constellation-filtered LOS satellites
filtered_ids = all_svid_prefix(ismember(all_constellations, sim_params_constellations));

filtered_los_sat_params = los_sats_params(contains(los_sats_params.Source, filtered_ids), :);
time_start.Format = 'HH:mm:ss';
time_end.Format = 'HH:mm:ss';

if any(filtered_los_sat_params.Duration < seconds(time_end - time_start))
    discarted_sats = filtered_los_sat_params(filtered_los_sat_params.Duration < seconds(time_end - time_start), :);
    for i = height(discarted_sats)
        log.info('The satellite %s was in LOS with the receiver, but its observation window of did not include the whole interval [%s, %s] and was therefore discarded.', discarted_sats(i,:).Source, time_start, time_end);
    end
    filtered_los_sat_params = filtered_los_sat_params(filtered_los_sat_params.Duration == seconds(time_end - time_start), :);
end

%% Get SVID-filtered LOS satellites
% TODO: also filter by SVIDs

end

