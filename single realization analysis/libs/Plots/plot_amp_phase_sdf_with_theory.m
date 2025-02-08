function plot_amp_phase_sdf_with_theory( ...
    irr_params, ...
    scint_field, ...
    detrended_phase_realization, ...
    norm_phase_sdf, ...
    doppler_frequency, ...
    mu)
% plot_amp_phase_sdf_with_theory
%
% Syntax:
%   plot_amp_phase_sdf_with_theory(irr_params, scint_field, ...
%                                  detrended_phase_realization, ...
%                                  norm_phase_sdf, ...
%                                  doppler_frequency, mu)
%
% Description:
%   Plots two subplots comparing:
%     1) The post-propagation intensity SDF from scint_field vs. the
%        theoretical intensity spectrum from Ispectrum.m.
%     2) The phase SDFs for:
%        - The pre-propagation phase (detrended_phase_realization),
%        - The post-propagation phase (phase(scint_field)),
%        - The theoretical phase SDF (norm_phase_sdf).
%
%   The code takes the one-sided (positive frequency) portion of the FFT, 
%   converts power to dB, and overlays the theoretical curves by interpolating 
%   them to the same frequency axis.
%
% Inputs:
%   irr_params               - Struct with fields .U, .p1, .p2, .mu0
%   scint_field              - 1D complex array (post-propagation)
%   detrended_phase_realization - 1D real array (pre-propagation phase)
%   norm_phase_sdf           - 1D array of theoretical phase SDF
%   doppler_frequency        - 1D array of length = nfft
%   mu                       - 1D array of length = nfft
%
% Outputs:
%   None (generates a figure with two subplots).
%
% Example:
%   irr_params = struct('U',1.5,'p1',2.4,'p2',3.7,'mu0',0.55);
%   plot_amp_phase_sdf_with_theory(irr_params, scint_field, ...
%                                  detrended_phase_realization, ...
%                                  norm_phase_sdf, ...
%                                  doppler_frequency, mu);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % (Implementation code as previously shown)
    % ----------------------------------------

    nfft = length(doppler_frequency);

    % Indices for positive frequencies (excluding DC)
    freq_idx_1sided = (nfft/2 + 2) : nfft;
    sdf_idx_1sided  = 2 : (nfft/2);
    partial_doppler_freq = doppler_frequency(freq_idx_1sided);
    partial_mu           = mu(freq_idx_1sided);

    % Helper for one-sided SDF
    compute_sdf_1sided = @(x) local_compute_sdf_1sided(x, nfft, doppler_frequency, sdf_idx_1sided);

    % Post-propagation intensity
    intensity_sdf_1sided_post = compute_sdf_1sided(abs(scint_field).^2);

    % Theoretical intensity from Ispectrum
    [Imu, muAxis, S4_val] = Ispectrum(irr_params.U, irr_params.p1, irr_params.p2, irr_params.mu0);
    logImu_interp = interp1(log10(muAxis), log10(Imu), log10(partial_mu), 'pchip', 'extrap');
    Imu_interp    = 10.^logImu_interp;

    % Phase SDF: pre- and post-propagation
    preprop_phase_sdf_1sided  = compute_sdf_1sided(detrended_phase_realization);
    postprop_phase_sdf_1sided = compute_sdf_1sided(phase(scint_field));

    % Shift the theoretical phase SDF
    theory_shifted = fftshift(norm_phase_sdf);
    phase_sdf_1sided_theory = theory_shifted(sdf_idx_1sided);

    % Plot
    figure('Name','Intensity & Phase SDF','Color',[1,1,1]);

    % Subplot 1: Intensity
    subplot(2,1,1); hold on; grid on
    plot(partial_doppler_freq, 10*log10(intensity_sdf_1sided_post), 'b', 'LineWidth',1.2)
    plot(partial_doppler_freq, 10*log10(Imu_interp), 'r--', 'LineWidth',1.2)
    set(gca,'XScale','log')
    xlabel('Doppler Frequency (Hz)')
    ylabel('SDF (dB)')
    title(sprintf('Intensity SDF (S4=%.3f)', S4_val))
    legend('Post-Prop','Theory','Location','best')

    % Subplot 2: Phase
    subplot(2,1,2); hold on; grid on
    plot(partial_doppler_freq, 10*log10(preprop_phase_sdf_1sided),  'b',  'LineWidth',1.2)
    plot(partial_doppler_freq, 10*log10(postprop_phase_sdf_1sided), 'g',  'LineWidth',1.2)
    plot(partial_doppler_freq, 10*log10(phase_sdf_1sided_theory),   'r--','LineWidth',1.2)
    set(gca,'XScale','log')
    xlabel('Doppler Frequency (Hz)')
    ylabel('SDF (dB)')
    title('Phase SDF')
    legend('Pre-Prop','Post-Prop','Theory','Location','best')

    sgtitle('Intensity & Phase SDF Comparison')

end

function sdf_1sided = local_compute_sdf_1sided(real_signal, nfft, doppler_freq, sdf_idx_1sided)
    raw_fft = abs(fft(real_signal, nfft)).^2 / nfft;
    sdf_partial = raw_fft(sdf_idx_1sided);
    df = abs(doppler_freq(2) - doppler_freq(1));
    sdf_1sided = sdf_partial / (nfft * df);
end
