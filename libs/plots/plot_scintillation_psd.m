function plot_scintillation_psd(cspsm_root_dir, out)
% plot_all_amp_phase_psds
%
% Syntax:
%   plot_all_amp_phase_psds(scint_field_struct, irr_params, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector)
%
% Description:
%   Generates **2x3 subplots** for each scintillation intensity level (Severe, Moderate, Weak).
%   - **Columns**: L1, L2, L5 frequency bands.
%   - **Rows**: Intensity PSD (top) and Phase PSD (bottom).
%
%   The intensity PSD of the post-propagation scintillation field is
%   compared with the theoretical intensity spectrum from `Ispectrum.m`.
%
%   The phase PSDs compare:
%   - Pre-propagation phase realization (`detrended_phase_realization`)
%   - Post-propagation phase (`phase(scint_field)`)
%   - Theoretical phase PSD (`norm_phase_psd`)
%
% Inputs:
%   scint_field_struct - Struct containing scintillation realizations for each scenario:
%       .Severe.L1, .Severe.L2, .Severe.L5
%       .Moderate.L1, .Moderate.L2, .Moderate.L5
%       .Weak.L1, .Weak.L2, .Weak.L5
%
%   irr_params - Struct with spectral parameters (.U, .p1, .p2, .mu0) for each scenario.
%   doppler_frequency_struct - Struct with Doppler frequency arrays for each scenario & freq.
%   mu_struct - Struct with the normalized wavenumber arrays for each scenario & freq.
%   rhof_veff_ratio_vector - 1x3 array containing (rho_F / v_eff) ratios for L1, L2, L5.
%
%
% Outputs:
%   None (generates multiple figures, one for each intensity level).
%
% Notes:
%   - This code is an adaptation of [1] and [2]
%
% References:
%   [1] "GenerateGPSPhaseScreenRealization.m" from the GNSS Scintillation
%       Simulator examples, available at
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator/blob/master/examples/GenerateGPSPhaseScreenRealization.m
%       [Accessed: 10-02-2025].
%
%   [2] "Display_SpectraModel.m" from the GNSS Scintillation Simulator
%       examples, available at
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator/blob/master/examples/Display_SpectraModel.m
%       [Accessed: 10-02-2025].
%
% Example:
%   plot_all_amp_phase_psds(scint_field_struct, irr_params_set, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

%% Initalize
severity = out.severity;

%% PSD of amplitude and phase of the ionospheric scintllation
% all constellation names
doppler_frequency = out.doppler_frequency_support;
constellations = setxor(string(fieldnames(out)), ...
    ["doppler_frequency_support", "time_utc", "satelliteScenario", "severity"]);
for constellation = constellations
    % Frequency names for this constellation
    freq_names = string(fieldnames(out.(constellation).spectral)).';
    % for all rx-sat scenario
    for i = 1:numel(out.(constellation).scenario)
        % Create new figure for each scenario
        fig = figure('Name', sprintf('%s Scintillation - Intensity & Phase PSDs', severity), 'Position',[50,50,1400,550]);
        set(fig, 'DefaultTextFontName', 'Helvetica');
        tiledlayout(2, numel(freq_names), "TileSpacing", "compact");
        sgtitle(sprintf('Scintillation PSD Analysis: Intensity & Phase (%s) for %s satellite %s', ...
            severity, upper(out.(constellation).scenario(i).sat.OrbitPropagator), ...
            out.(constellation).scenario(i).sat.Name), 'FontSize', 30, ...
            'FontName', 'Helvetica');
        for j = 1:numel(freq_names)
            % get parameters
            freq_name = freq_names(j);
            spectral_params = out.(constellation).spectral.(freq_name);
            mu = out.(constellation).scenario(i).(freq_name).mu;
            norm_phase_psd = out.(constellation).scenario(i).(freq_name).norm_phase_psd;
            detrended_phase = out.(constellation).scenario(i).(freq_name).detrended_phase;
            scint_field = out.(constellation).scenario(i).(freq_name).scint_field;
            rhof_veff_ratio = out.(constellation).scenario(i).(freq_name).rhof_veff_ratio;
            % Compute one-sided Doppler frequencies
            nfft = length(doppler_frequency);
            freq_idx_1sided = (nfft/2 + 2) : nfft;
            psd_idx_1sided  = 2 : (nfft/2);
            partial_doppler_freq = doppler_frequency(freq_idx_1sided);
            partial_mu           = mu(freq_idx_1sided);
            
            % Compute Simulated Intensity PSD (Post-Propagation)
            intensity_signal = abs(scint_field).^2;
            intensity_psd_1sided_post = compute_psd_1sided(intensity_signal, nfft, doppler_frequency, psd_idx_1sided);
            experimental_s4_val = get_S4(intensity_signal);

            % Compute Simulated Post- and Pre-Propagated Phase PSDs
            % NOTE: Both intensity and phase are normalized by
            % `rhof_veff_ratio` to agree with the code in [2].

            % NOTE: `detrended_phase` contains only the refractive-related
            % effect of the phase disturbance at the IPP point,
            % which has not been propagated to the receiver yet
            preprop_phase_psd_1sided  = compute_psd_1sided(detrended_phase, nfft, doppler_frequency, psd_idx_1sided) / rhof_veff_ratio;
            % NOTE: `phase(scint_field)` is the phase of the complex field
            % after the propagation, which contains noy only the refractive
            % part, but also the difracted part caused by the free-space
            % propagation
            postprop_phase_psd_1sided = compute_psd_1sided(phase(scint_field), nfft, doppler_frequency, psd_idx_1sided) / rhof_veff_ratio;

            % Compute Theoretical Intensity PSD (Interpolated)
            [Imu, muAxis, theoretical_s4_val] = Ispectrum(cspsm_root_dir, ...
                spectral_params.U, spectral_params.p1, spectral_params.p2, ...
                spectral_params.mu0);
            logImu_interp = interp1(log10(muAxis), log10(Imu), log10(partial_mu), 'pchip', 'extrap');
            Imu_interp = 10.^logImu_interp;
            % Compute Theoretical phase PSD
            phase_psd_1sided_theory = norm_phase_psd(freq_idx_1sided);

            % Plot Simulated Intensity PSD
            nexttile(j);
            hold on;
            plot(partial_doppler_freq, 10*log10(intensity_psd_1sided_post), 'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(Imu_interp), 'r--', 'LineWidth', 1.2);
            set(gca, 'XScale', 'log', 'FontName', 'Helvetica', 'FontSize', 20);
            xlabel('Doppler Frequency (Hz)', 'FontName', 'Helvetica', 'FontSize', 10);
            ylabel('Intensity PSD [dB]', 'FontName', 'Helvetica', 'FontSize', 10);
            title(sprintf('%s %s | S_4: Th. %.3f, Exp. %.3f | ρ_F/v_{eff} = %.2f', severity, freq_name, theoretical_s4_val, experimental_s4_val, rhof_veff_ratio), 'FontSize', 20, 'FontName', 'Helvetica');
            legend( ...
                'Post-Propagation PSD', ...
                'Theoretical Model', ...
                'Location', 'best', ...
                'FontName', 'Helvetica', ...
                'FontSize', 15);
            grid on;
            grid minor;
            hold off;

            % Plot Phase PSD
            nexttile(j + numel(freq_names));
            hold on;
            plot(partial_doppler_freq, 10*log10(preprop_phase_psd_1sided), ...
                'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(postprop_phase_psd_1sided), ...
                'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(phase_psd_1sided_theory), ...
                'LineWidth', 1.2);
            set(gca, 'XScale', 'log', 'FontName', 'Helvetica');
            xlabel('Doppler Frequency (Hz)', 'FontName', 'Helvetica');
            ylabel('Phase PSD [dB]', 'FontName', 'Helvetica');
            title(sprintf('%s %s Phase PSD | ρ_F/v_{eff} = %.2f', severity, freq_name, rhof_veff_ratio), ...
                'FontSize', 20, 'FontName', 'Helvetica');
            legend('Pre-Propagation PSD', 'Post-Propagation PSD', ...
                'Theoretical Model', 'Location', 'best', 'FontName', ...
                'Helvetica', 'FontSize', 15);
            grid on;
            grid minor;
            hold off;
        end
    end
end



% -------------------------------------------------------------------------

    function psd_1sided = compute_psd_1sided(real_signal, nfft, doppler_freq, psd_idx_1sided)
        % Compute the one-sided power spectral density function (PSD)
        raw_fft = abs(fft(real_signal, nfft)).^2 / nfft;
        partial_psd = raw_fft(psd_idx_1sided);
        df = abs(doppler_freq(2) - doppler_freq(1));
        psd_1sided = partial_psd / (nfft * df);
    end

end