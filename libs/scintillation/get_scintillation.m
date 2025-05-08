function out = get_scintillation(sim_params)

%% Initialize
rx = sim_params.satelliteScenario.Platforms(1);

%% Compute the nongeometric-dependent outputs
% They are
% - U and μ₀: dependent only on the frequency (computed via extrapolation)
% - p₁ and p₂: independent
% - Doppler frequency: dependent only on t_samp and FFT samples
% Frequency support where the PSD are plotted
nfft = nicefftnum(sim_params.sim_time / sim_params.t_samp);
out.doppler_frequency_support = (-nfft/2 : nfft/2-1) / (nfft * sim_params.t_samp);
for constellation = sim_params.constellations
    % initialize a 0x0 struct array which will store geometric-related
    % outputs. Each, struct scalar in this array is related to a sat-rx
    % geometry
    out.(constellation).scenario = struct([]);
    for i = 1:numel(sim_params.freqs.(constellation).value)
        freq_name = sim_params.freqs.(constellation).name(i);
        freq_value = sim_params.freqs.(constellation).value(i);
        
        % extrapolate U, and μ₀ from `freq_ref` to `freq`
        out.(constellation).spectral.(freq_name) = ...
            extrapolate_spectral_params(sim_params, freq_value);
    end
end

%% Compute the geometric related output
% for all filtered LOS sat (of a given contellation)
for sat = sim_params.satelliteScenario.Satellites
    % get LOS sat and constellation name
    sat_constellation = sat.OrbitPropagator;
    % save sat and receiver to output
    out.(sat_constellation).scenario(end+1).sat = sat;
    out.(sat_constellation).scenario(end).rx = rx;
    % for this geometry propagation, get the ρF/veff for the reference
    % frequency
    rhof_veff_ratio_ref = get_scaling_param(rx, sat, sim_params);

    % for all valid frequencies for this satellite contellation
    for j = 1:numel(sim_params.freqs.(sat_constellation).name)
        freq_name = sim_params.freqs.(sat_constellation).name(j);
        freq_value = sim_params.freqs.(sat_constellation).value(j);

        rhof_veff_ratio = extrapolate_scaling_param(sim_params, ...
            freq_value, rhof_veff_ratio_ref);
        
        % get scintillation realization
        [scint_field, theo_phase_psd, detrended_phase, mu, ...
            postprop_amplitude, postprop_phase, ...
            intensity_psd_1sided_post, s4, ...
            preprop_phase_psd_1sided, postprop_phase_psd_1sided] = ...
            get_scintillation_realization( ...
            sim_params, ...
            out, ...
            constellation, ...
            freq_name, ...
            rhof_veff_ratio, ...
            nfft);
        
        %% outputs
        % frequency support
        out.(sat_constellation).scenario(end).(freq_name).mu = mu;
        out.(sat_constellation).scenario(end).(freq_name).complex_field_postprop = scint_field;
        % Computed S4
        out.(sat_constellation).scenario(end).(freq_name).S4 = s4;
        % amplitude
        out.(sat_constellation).scenario(end).(freq_name).amplitude.timeseries_postprop = postprop_amplitude;
        out.(sat_constellation).scenario(end).(freq_name).amplitude.psd_postprop = intensity_psd_1sided_post;
        % phase
        out.(sat_constellation).scenario(end).(freq_name).phase.timeseries.postprop = postprop_phase;
        out.(sat_constellation).scenario(end).(freq_name).phase.timeseries.preprop = detrended_phase;
        out.(sat_constellation).scenario(end).(freq_name).phase.psd.theo_phase = theo_phase_psd;
        out.(sat_constellation).scenario(end).(freq_name).phase.psd.preprop_1sided = preprop_phase_psd_1sided;
        out.(sat_constellation).scenario(end).(freq_name).phase.psd.postprop_1sided = postprop_phase_psd_1sided;
        % ρF/veff
        out.(sat_constellation).scenario(end).(freq_name).rhof_veff_ratio = rhof_veff_ratio;
    end
end

%% Get output directly from simulation parameters
% Pass `satelliteScenario` object to output
out.satelliteScenario = sim_params.satelliteScenario;
% Severity
out.severity = sim_params.severity;
end

