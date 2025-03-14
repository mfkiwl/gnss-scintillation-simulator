%% Evaluate Precision Loss in Complex Scintillation Time Series (L1, Severe Case)
% This script generates a complex scintillation time series for the L1
% frequency (Severe case) using your scintillation model, then converts the
% series from double precision to single and half precision. It computes error
% metrics for both amplitude and phase and produces histograms and sample
% comparisons for both quantities.
%
% Author:
%   Your Name
%   Email: your.email@example.com
%   Date: YYYY-MM-DD

clearvars; clc;

%% Add required paths (adjust as needed)
addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Model Setup & Time Series Generation (Severe, L1)
general_params = get_general_parameters();
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);
irr_params_set = get_irregularity_parameters();
[extrapolated_irr_params.Severe, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params_set.Severe, general_params, rhof_veff_ratio_L1);

seed = 1;
[scint_L1, ~, ~, ~, ~] = get_scintillation_time_series( ...
    general_params, ...
    extrapolated_irr_params.Severe.L1, ...
    rhof_veff_ratio_vector(1), ...
    seed);
% scint_L1 is a complex double time series

%% Convert to Different Precision Formats
scint_L1_single = single(scint_L1);
scint_L1_half   = half(scint_L1);  % Requires MATLAB R2023a+

%% Compute Error Metrics for Amplitude and Phase
% Original amplitude and phase
amp_orig   = abs(scint_L1);
phase_orig = angle(scint_L1);

% Single precision errors
amp_single   = abs(scint_L1_single);
phase_single = angle(scint_L1_single);
amp_error_single   = amp_orig - amp_single;
phase_error_single = angle(exp(1i*(phase_orig - phase_single)));  % wrapped phase error

% Half precision errors
amp_half   = abs(scint_L1_half);
phase_half = angle(scint_L1_half);
amp_error_half   = amp_orig - amp_half;
phase_error_half = angle(exp(1i*(phase_orig - phase_half)));

% Compute error metrics for amplitude
mse_amp_single    = mean(abs(amp_error_single).^2);
max_amp_error_single = max(abs(amp_error_single));

mse_amp_half    = mean(abs(amp_error_half).^2);
max_amp_error_half = max(abs(amp_error_half));

% Compute error metrics for phase
mse_phase_single    = mean(abs(phase_error_single).^2);
max_phase_error_single = max(abs(phase_error_single));

mse_phase_half    = mean(abs(phase_error_half).^2);
max_phase_error_half = max(abs(phase_error_half));

%% Display Error Metrics
fprintf('=== Amplitude Error Metrics ===\n');
fprintf('Single Precision: MSE = %e, Max Error = %e\n', mse_amp_single, max_amp_error_single);
fprintf('Half Precision  : MSE = %e, Max Error = %e\n\n', mse_amp_half, max_amp_error_half);

fprintf('=== Phase Error Metrics ===\n');
fprintf('Single Precision: MSE = %e, Max Error = %e\n', mse_phase_single, max_phase_error_single);
fprintf('Half Precision  : MSE = %e, Max Error = %e\n\n', mse_phase_half, max_phase_error_half);

%% Plot Histograms of Error Magnitudes
figure;
subplot(2,2,1);
histogram(double(abs(amp_error_single)), 100);
title('Amplitude Error: Single');
xlabel('Error'); ylabel('Frequency');

subplot(2,2,2);
histogram(double(abs(amp_error_half)), 100);
title('Amplitude Error: Half');
xlabel('Error'); ylabel('Frequency');

subplot(2,2,3);
histogram(double(abs(phase_error_single)), 100);
title('Phase Error: Single');
xlabel('Error (rad)'); ylabel('Frequency');

subplot(2,2,4); 
histogram(double(abs(phase_error_half)), 100);
title('Phase Error: Half');
xlabel('Error (rad)'); ylabel('Frequency');

%% Plot Sample Comparisons of Amplitude and Phase
% Define a sample range for plotting (e.g., samples 10000 to 10100)
sample_range = 10000:10100;

figure;
tiledlayout(2,1, 'TileSpacing','compact','Padding','compact');

% Plot Amplitude Comparison
nexttile;
plot(sample_range, abs(scint_L1(sample_range)), 'b-', 'LineWidth', 1.5); hold on;
plot(sample_range, abs(scint_L1_single(sample_range)), 'r--', 'LineWidth', 1.5);
plot(sample_range, abs(scint_L1_half(sample_range)), 'g-.', 'LineWidth', 1.5);
title('Amplitude Comparison');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Original (Double)', 'Single', 'Half', 'Location','best');
hold off;

% Plot Phase Comparison
nexttile;
plot(sample_range, angle(scint_L1(sample_range)), 'b-', 'LineWidth', 1.5); hold on;
plot(sample_range, angle(scint_L1_single(sample_range)), 'r--', 'LineWidth', 1.5);
plot(sample_range, angle(scint_L1_half(sample_range)), 'g-.', 'LineWidth', 1.5);
title('Phase Comparison');
xlabel('Sample Index');
ylabel('Phase (rad)');
legend('Original (Double)', 'Single', 'Half', 'Location','best');
hold off;
