function filtered_los_sat_params = get_filtered_los_sat_params(log, sim_params, los_sats_params)
%GET_LOS_SATS Summary of this function goes here
%   Detailed explanation goes here
%% Initialization
all_constellations = sim_params.cte.all_constellations;
% sim_params_svids = sim_params.svids;
sim_params_constellations = sim_params.constellations;
all_svid_prefix = sim_params.cte.all_svid_prefix;
trange = sim_params.time_range;

%% Get constellation-filtered LOS satellites
filtered_ids = all_svid_prefix(contains(all_constellations, sim_params_constellations));

filtered_los_sat_params = los_sats_params(contains(los_sats_params.Source, filtered_ids), :);
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
% TODO: also filter by SVIDs

end

