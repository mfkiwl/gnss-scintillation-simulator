function plot_scintillation_realization(scint_realization, severity, cspsm_root_dir)
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

%% PSD of amplitude and phase of the ionospheric scintllation

% all constellation names
constellations = string(fieldnames(scint_realization));
for constellation = constellations
    % Frequency names for this constellation
    freq_names = setxor(string(fieldnames(scint_realization.(constellation))), ["sat", "rx"]).';
    % Create new figure for each scenario
    fig = figure('Name', sprintf('%s Scintillation - Intensity & Phase PSDs', severity), 'Position',[50,50,1400,550]);
    set(fig, 'DefaultTextFontName', 'Helvetica');
    tiledlayout(2, numel(freq_names), "TileSpacing", "compact");
    for j = 1:numel(freq_names)
        freq_name = freq_names(j);
        % get scintillation realization parameters for this sat-rx scenario
        % (get only for the very first)
        scint_field = scint_realization.(constellation)(1).(freq_name).scint_field;
        detrended_phase = scint_realization.(constellation)(1).(freq_name).detrended_phase;
        spectral_params = scint_realization.(constellation)(1).(freq_name).spectral_params;
        rhof_veff_ratio = scint_realization.(constellation)(1).(freq_name).rhof_veff_ratio;
        doppler_frequency = scint_realization.(constellation)(1).(freq_name).doppler_frequency;
        mu = scint_realization.(constellation)(1).(freq_name).mu;

        % Compute one-sided Doppler frequencies
        nfft = length(doppler_frequency);
        freq_idx_1sided = (nfft/2 + 2) : nfft;
        psd_idx_1sided  = 2 : (nfft/2);
        partial_doppler_freq = doppler_frequency(freq_idx_1sided);
        partial_mu           = mu(freq_idx_1sided);

        % Compute Intensity PSD (Post-Propagation)
        intensity_signal = abs(scint_field).^2;
        intensity_psd_1sided_post = compute_psd_1sided(intensity_signal, nfft, doppler_frequency, psd_idx_1sided);
        experimental_s4_val = get_S4(intensity_signal);

        % Compute Theoretical Intensity PSD (Interpolated)
        [Imu, muAxis, theoretical_s4_val] = Ispectrum(cspsm_root_dir, spectral_params.U, spectral_params.p1, spectral_params.p2, spectral_params.mu0);
        logImu_interp = interp1(log10(muAxis), log10(Imu), log10(partial_mu), 'pchip', 'extrap');
        Imu_interp = 10.^logImu_interp;

        % Compute Phase PSDs
        % NOTE: Both intensity and phase are normalized by
        % `rhof_veff_ratio` to agree with the code in [2].
        preprop_phase_psd_1sided  = compute_psd_1sided(detrended_phase, nfft, doppler_frequency, psd_idx_1sided) / rhof_veff_ratio;
        postprop_phase_psd_1sided = compute_psd_1sided(phase(scint_field), nfft, doppler_frequency, psd_idx_1sided) / rhof_veff_ratio;
        norm_phase_psd = get_norm_phase_psd(mu, spectral_params);
        phase_psd_1sided_theory = norm_phase_psd(freq_idx_1sided);

        % Plot Intensity PSD
        nexttile(j);
        hold on;
        plot(partial_doppler_freq, 10*log10(intensity_psd_1sided_post), 'b', 'LineWidth', 1.2);
        plot(partial_doppler_freq, 10*log10(Imu_interp), 'r--', 'LineWidth', 1.2);
        set(gca, 'XScale', 'log', 'FontName', 'Helvetica');
        xlabel('Doppler Frequency (Hz)', 'FontName', 'Helvetica');
        ylabel('Intensity PSD [dB]', 'FontName', 'Helvetica');
        title(sprintf('%s %s | S_4: Th. %.3f, Exp. %.3f', severity, freq_name, theoretical_s4_val, experimental_s4_val), 'FontSize', 10, 'FontName', 'Helvetica');
        legend('Post-Propagation PSD', 'Theoretical Model', 'Location', 'best', 'FontName', 'Helvetica');
        grid on;
        grid minor;
        hold off;

        % Plot Phase PSD
        nexttile(j + numel(freq_names));
        hold on;
        plot(partial_doppler_freq, 10*log10(preprop_phase_psd_1sided), 'b', 'LineWidth', 1.2);
        plot(partial_doppler_freq, 10*log10(postprop_phase_psd_1sided), 'g', 'LineWidth', 1.2);
        plot(partial_doppler_freq, 10*log10(phase_psd_1sided_theory), 'r--', 'LineWidth', 1.2);
        set(gca, 'XScale', 'log', 'FontName', 'Helvetica');
        xlabel('Doppler Frequency (Hz)', 'FontName', 'Helvetica');
        ylabel('Phase PSD [dB]', 'FontName', 'Helvetica');
        title(sprintf('%s %s Phase PSD', severity, freq_name), 'FontSize', 10, 'FontName', 'Helvetica');
        legend('Pre-Propagation PSD', 'Post-Propagation PSD', 'Theoretical Model', 'Location', 'best', 'FontName', 'Helvetica');
        grid on;
        grid minor;
        hold off;
    end
    sgtitle(sprintf('Scintillation PSD Analysis: Intensity & Phase (%s) for %s satellite %s', severity, upper(scint_realization.(constellation)(1).sat.OrbitPropagator), scint_realization.(constellation)(1).sat.Name), 'FontSize', 12, 'FontName', 'Helvetica');
end

    function psd_1sided = compute_psd_1sided(real_signal, nfft, doppler_freq, psd_idx_1sided)
        % Compute the one-sided power spectral density function (PSD)
        raw_fft = abs(fft(real_signal, nfft)).^2 / nfft;
        partial_psd = raw_fft(psd_idx_1sided);
        df = abs(doppler_freq(2) - doppler_freq(1));
        psd_1sided = partial_psd / (nfft * df);
    end

end