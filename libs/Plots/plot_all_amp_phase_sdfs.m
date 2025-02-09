function plot_all_amp_phase_sdfs(scint_field_struct, irr_params, detrended_phase_realization_struct, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector)
% plot_all_amp_phase_sdfs
%
% Syntax:
%   plot_all_amp_phase_sdfs(scint_field_struct, irr_params, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector)
%
% Description:
%   Generates **2x3 subplots** for each scintillation intensity level (Severe, Moderate, Mild).
%   - **Columns**: L1, L2, L5 frequency bands.
%   - **Rows**: Intensity SDF (top) and Phase SDF (bottom).
%
%   The intensity SDF of the **post-propagation** scintillation field is compared with 
%   the theoretical intensity spectrum from `Ispectrum.m`.
%   
%   The phase SDFs compare:
%   - Pre-propagation phase realization (`detrended_phase_realization`)
%   - Post-propagation phase (`phase(scint_field)`)
%   - Theoretical phase SDF (`norm_phase_sdf`)
%
% Inputs:
%   scint_field_struct - Struct containing scintillation realizations for each scenario:
%       .Severe.L1, .Severe.L2, .Severe.L5
%       .Moderate.L1, .Moderate.L2, .Moderate.L5
%       .Mild.L1, .Mild.L2, .Mild.L5
%
%   irr_params - Struct with spectral parameters (.U, .p1, .p2, .mu0) for each scenario.
%   doppler_frequency_struct - Struct with Doppler frequency arrays for each scenario & freq.
%   mu_struct - Struct with the normalized wavenumber arrays for each scenario & freq.
%   rhof_veff_ratio_vector - 1x3 array containing (rho_F / v_eff) ratios for L1, L2, L5.
%
% Outputs:
%   None (generates multiple figures, one for each intensity level).
%
% Example:
%   plot_all_amp_phase_sdfs(scint_field_struct, irr_params_set, doppler_frequency_struct, mu_struct, rhof_veff_ratio_vector);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    scenarios = fieldnames(scint_field_struct);  % {'Severe', 'Moderate', 'Mild'}
    frequencies = {'L1', 'L2', 'L5'};

    for i = 1:numel(scenarios)
        scenario = scenarios{i};
        rhof_veff_ratio = rhof_veff_ratio_vector(i);
        
        % Create new figure for each scenario
        fig = figure('Name', sprintf('%s Scintillation - Intensity & Phase SDFs', scenario), 'Position',[50,50,1200,700]);
        set(fig, 'DefaultTextFontName', 'Helvetica');
        tiledlayout(2, 3, "TileSpacing", "compact");

        for j = 1:numel(frequencies)
            freq = frequencies{j};
            scint_field = scint_field_struct.(scenario).(freq);
            detrended_phase_realization = detrended_phase_realization_struct.(scenario).(freq);
            irr_params_freq = irr_params.(scenario).(freq);
            doppler_frequency = doppler_frequency_struct.(scenario).(freq);
            mu = mu_struct.(scenario).(freq);

            % Compute one-sided Doppler frequencies
            nfft = length(doppler_frequency);
            freq_idx_1sided = (nfft/2 + 2) : nfft; 
            sdf_idx_1sided  = 2 : (nfft/2);        
            partial_doppler_freq = doppler_frequency(freq_idx_1sided);
            partial_mu           = mu(freq_idx_1sided);

            % Compute Intensity SDF (Post-Propagation)
            intensity_signal = abs(scint_field).^2;
            intensity_sdf_1sided_post = compute_sdf_1sided(intensity_signal, nfft, doppler_frequency, sdf_idx_1sided);
            experimental_s4_val = get_S4(intensity_signal);

            % Compute Theoretical Intensity SDF (Interpolated)
            [Imu, muAxis, theoretical_s4_val] = Ispectrum(irr_params_freq.U, irr_params_freq.p1, irr_params_freq.p2, irr_params_freq.mu0);
            logImu_interp = interp1(log10(muAxis), log10(Imu), log10(partial_mu), 'pchip', 'extrap');
            Imu_interp = 10.^logImu_interp;

            % Compute Phase SDFs
            preprop_phase_sdf_1sided  = compute_sdf_1sided(detrended_phase_realization, nfft, doppler_frequency, sdf_idx_1sided) / rhof_veff_ratio;
            postprop_phase_sdf_1sided = compute_sdf_1sided(phase(scint_field), nfft, doppler_frequency, sdf_idx_1sided) / rhof_veff_ratio;
            norm_phase_sdf = get_norm_phase_sdf(mu, irr_params_freq);
            phase_sdf_1sided_theory = norm_phase_sdf(freq_idx_1sided);

            % Plot Intensity SDF
            nexttile(j);
            hold on;
            plot(partial_doppler_freq, 10*log10(intensity_sdf_1sided_post), 'b', 'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(Imu_interp), 'r--', 'LineWidth', 1.2);
            set(gca, 'XScale', 'log', 'FontName', 'Helvetica');
            xlabel('Frequency (Hz)', 'FontName', 'Helvetica');
            ylabel('Intensity SDF [dB]', 'FontName', 'Helvetica');
            title(sprintf('%s %s | S4: Th. %.3f, Exp. %.3f', scenario, freq, theoretical_s4_val, experimental_s4_val), 'FontSize', 10, 'FontName', 'Helvetica');
            legend('Post-Propagation SDF', 'Theoretical Model', 'Location', 'best', 'FontName', 'Helvetica');
            grid on;
            grid minor;
            hold off;

            % Plot Phase SDF
            nexttile(j + 3);
            hold on;
            plot(partial_doppler_freq, 10*log10(preprop_phase_sdf_1sided), 'b', 'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(postprop_phase_sdf_1sided), 'g', 'LineWidth', 1.2);
            plot(partial_doppler_freq, 10*log10(phase_sdf_1sided_theory), 'r--', 'LineWidth', 1.2);
            set(gca, 'XScale', 'log', 'FontName', 'Helvetica');
            xlabel('Frequency (Hz)', 'FontName', 'Helvetica');
            ylabel('Phase SDF [dB]', 'FontName', 'Helvetica');
            title(sprintf('%s %s Phase SDF', scenario, freq), 'FontSize', 10, 'FontName', 'Helvetica');
            legend('Pre-Propagation SDF', 'Post-Propagation SDF', 'Theoretical Model', 'Location', 'best', 'FontName', 'Helvetica');
            grid on;
            grid minor;
            hold off;
        end

        sgtitle(sprintf('Scintillation SDF Analysis: Intensity & Phase (%s)', scenario), 'FontSize', 12, 'FontName', 'Helvetica');

        % Copy figure size to paper size for correct PDF export in centimeters
        set(fig, 'PaperUnits', 'centimeters');
        fig_pos = get(fig, 'Position'); % Get figure size in pixels
        fig_width_cm = fig_pos(3) * 2.54 / 96;  % Convert from pixels to cm
        fig_height_cm = fig_pos(4) * 2.54 / 96;
        set(fig, 'PaperSize', [fig_width_cm fig_height_cm]);  
        set(fig, 'PaperPosition', [0 0 fig_width_cm fig_height_cm]); % Match figure size exactly

        % Save as PDF (high-quality, uncompressed)
        pdf_filename = sprintf('%s_Scintillation_SDFs.pdf', scenario);
        print(fig, fullfile('..','figures',pdf_filename), '-dpdf', '-vector');
        
        disp(['Saved: ', fullfile('..','figures',pdf_filename)]);
    end
end

function sdf_1sided = compute_sdf_1sided(real_signal, nfft, doppler_freq, sdf_idx_1sided)
    % Compute the one-sided power spectral density function (SDF)
    raw_fft = abs(fft(real_signal, nfft)).^2 / nfft;
    partial_sdf = raw_fft(sdf_idx_1sided);
    df = abs(doppler_freq(2) - doppler_freq(1));
    sdf_1sided = partial_sdf / (nfft * df);
end
