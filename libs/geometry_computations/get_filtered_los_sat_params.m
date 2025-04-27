function filtered_los_sat_params = get_filtered_los_sat_params(log, all_constellations, sim_params_svids, sim_params_constellations, all_ids, trange, los_sats_params)
%GET_LOS_SATS Summary of this function goes here
%   Detailed explanation goes here

%% Get constellation-filtered LOS satellites
filtered_ids = all_ids(contains(all_constellations, sim_params_constellations));

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

