%% Main Simulation Script for GNSS Scintillation Dataset Generation
%
% Syntax:
%   run_main_simulation_script
%
% Description:
%   This script generates GNSS scintillation datasets by sweeping over various 
%   simulation parameters (drift velocities, CN0 values, severity levels, Monte Carlo
%   runs, and frequency bands). For each parameter combination, the script:
%       - Generates the scintillation time series.
%       - Computes the normalized autocorrelation of the scintillation intensity to 
%         classify the scintillation duration ('long', 'medium', or 'short').
%       - Adds thermal noise.
%       - Saves the received signal (split into real and imaginary parts) into an 
%         HDF5 file (organized by receiver city) with associated attributes.
%
%   This version prints high-level progress messages once per city. After each city,
%   it prints the elapsed time for that city, an estimate of the remaining time for the
%   simulation, and the overall progress percentage. A final summary log displays the
%   total elapsed time and average time per city.
%
% Example:
%   run_main_simulation_script
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

%% Initialization
clearvars; clc;
addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Parameter Setup
drift_velocities_amount = 1;  % number of drift velocities (e.g., from 25 to 125 m/s)
carrier_to_noise_ratios_amount = 1; % number of CN0 values (e.g., from 25 to 45 dB-Hz)
monte_carlo_runs = 1;
mask_angle_deg = 15;  % mask angle (degrees)

data_set_params = get_dataset_params(drift_velocities_amount, carrier_to_noise_ratios_amount, mask_angle_deg);
cities_str_array = string(fieldnames(data_set_params.specific));
severity_str_array = string(fieldnames(data_set_params.general.irr_params));

receiver_sampling_frequency = 2e7;
B = receiver_sampling_frequency / 2;
rx_mean_power = 1;
frequency_str_array = ["L1", "L2", "L5"];

filename = "simplified_cpsm_dataset.h5";

% Global HDF5 file attributes.
hdf5_file_attributes = struct(...
    'date_time', data_set_params.general.date_time, ...
    'speed_of_light', data_set_params.general.c, ...
    'gps_bands', data_set_params.general.gps_bands, ...
    'sampling_interval_after_correlation', data_set_params.general.dt, ...
    'simulation_time', data_set_params.general.simulation_time, ...
    'ionospheric_piercing_point_height', data_set_params.general.ipp_height, ...
    'drift_velocities_amount', data_set_params.general.drift_velocities_amount, ...
    'drift_velocities', data_set_params.general.drift_velocities(2,:), ...
    'carrier_to_noise_ratios_amount', data_set_params.general.carrier_to_noise_ratios_amount, ...
    'carrier_to_noise_ratios', data_set_params.general.carrier_to_noise_ratios);

% Initialize HDF5 file (with versioning/overwrite handling).
new_filename = initialize_hdf5_file(filename, cities_str_array.', hdf5_file_attributes);

% Seed counter for random number generation.
seed_counter = 0;

%% Compute Total Iterations (for overall progress)
total_iter = 0;
for city = cities_str_array.'
    prn_list = data_set_params.specific.(char(city)).prn_string_array;
    total_iter = total_iter + numel(prn_list) * drift_velocities_amount * numel(severity_str_array) * ...
                 numel(data_set_params.general.carrier_to_noise_ratios) * monte_carlo_runs;
end
overall_iter = 0;
overall_tic = tic;  % Start overall timer

fprintf('Starting simulation parameter sweeps...\n');

%% Main Loop: Process Each City
city_times = zeros(numel(cities_str_array),1);  % Store elapsed time per city
city_count = 0;
for city_str = cities_str_array.'
    city_count = city_count + 1;
    city_tic = tic;  % Start timer for current city
    
    % Set city-level attributes.
    city_attributes = struct( ...
        'rx_init_llh_radians', data_set_params.specific.(city_str).rx_init_llh, ...
        'rx_velocity_m_per_sec', data_set_params.specific.(city_str).rx_velocity, ...
        'prn_str_array', data_set_params.specific.(city_str).prn_string_array);
    set_city_group_attributes(new_filename, city_str, city_attributes);
    
    % Process each PRN for the city.
    prn_list = data_set_params.specific.(char(city_str)).prn_string_array;
    for prn_str = prn_list
        % (No per-PRN print to avoid overprinting.)
        for drift_idx = 1:drift_velocities_amount
            for severity_str = severity_str_array.'
                for carrier_to_noise_ratio = data_set_params.general.carrier_to_noise_ratios
                    for mc_run = 1:monte_carlo_runs
                        num_samples = data_set_params.general.simulation_time / data_set_params.general.dt;
                        real_received_signal = zeros(num_samples, numel(frequency_str_array), 'single');
                        imag_received_signal = zeros(num_samples, numel(frequency_str_array), 'single');
                        labels_struct = struct('label_L1', string(), 'label_L2', string(), 'label_L5', string());
                        
                        for frequency_str = frequency_str_array
                            frequency_idx = find(frequency_str_array == frequency_str);
                            
                            thermal_noise = get_thermal_noise( ...
                                data_set_params.general.simulation_time, ...
                                data_set_params.general.dt, ...
                                rx_mean_power, ...
                                carrier_to_noise_ratio, ...
                                B, ...
                                'data_type', 'single');
                            
                            U_str = ['U_', char(frequency_str)];
                            mu0_str = ['mu0_', char(frequency_str)];
                            irr_params.U = data_set_params.general.irr_params.(severity_str).(U_str);
                            irr_params.mu0 = data_set_params.general.irr_params.(severity_str).(mu0_str);
                            irr_params.p1 = data_set_params.general.irr_params.(severity_str).p1;
                            irr_params.p2 = data_set_params.general.irr_params.(severity_str).p2;
                            
                            rhof_veff_ratio_str = ['rhof_veff_ratio_', char(frequency_str)];
                            rhof_veff_ratio = data_set_params.specific.(city_str).rhof_veff_table{drift_idx, prn_str}.(rhof_veff_ratio_str);
                            
                            scint_series = get_scintillation_time_series( ...
                                data_set_params.general, ...
                                irr_params, ...
                                rhof_veff_ratio, ...
                                seed_counter, ...
                                'data_type', 'single');
                            scint_series_sliced = scint_series(1:num_samples).';
                            
                            % Compute label for current frequency.
                            label_field = ['label_', char(frequency_str)];
                            labels_struct.(label_field) = get_label_scint(scint_series_sliced, data_set_params.general.dt, ...
                                'lower_time', 0.5, 'upper_time', 2.0);
                            
                            received_signal = scint_series_sliced + thermal_noise;
                            
                            real_received_signal(:, frequency_idx) = real(received_signal);
                            imag_received_signal(:, frequency_idx) = imag(received_signal);
                        end
                        
                        % Build attributes for the current dataset.
                        dataset_attributes = struct( ...
                            'label_L1', labels_struct.label_L1, ...
                            'label_L2', labels_struct.label_L2, ...
                            'label_L5', labels_struct.label_L5, ...
                            'prn', prn_str, ...
                            'drift_velocity', data_set_params.general.drift_velocities(2, drift_idx), ...
                            'carrier_to_noise_ratio', carrier_to_noise_ratio, ...
                            'severity', severity_str, ...
                            'monte_carlo_run', mc_run, ...
                            'rhof_veff_ratio_L1', data_set_params.specific.(city_str).rhof_veff_table{drift_idx, prn_str}.rhof_veff_ratio_L1, ...
                            'rhof_veff_ratio_L2', data_set_params.specific.(city_str).rhof_veff_table{drift_idx, prn_str}.rhof_veff_ratio_L2, ...
                            'rhof_veff_ratio_L5', data_set_params.specific.(city_str).rhof_veff_table{drift_idx, prn_str}.rhof_veff_ratio_L5, ...
                            'p1', data_set_params.general.irr_params.(severity_str).p1, ...
                            'p2', data_set_params.general.irr_params.(severity_str).p2, ...
                            'U_L1', data_set_params.general.irr_params.(severity_str).U_L1, ...
                            'mu0_L1', data_set_params.general.irr_params.(severity_str).mu0_L1, ...
                            'U_L2', data_set_params.general.irr_params.(severity_str).U_L2, ...
                            'mu0_L2', data_set_params.general.irr_params.(severity_str).mu0_L2, ...
                            'U_L5', data_set_params.general.irr_params.(severity_str).U_L5, ...
                            'mu0_L5', data_set_params.general.irr_params.(severity_str).mu0_L5);
                        
                        % Write dataset to the HDF5 file.
                        write_received_signal_to_hdf5(new_filename, city_str, ...
                            real_received_signal, imag_received_signal, ...
                            dataset_attributes, 'chunk_size', [10000, 3], 'deflate_level', 4);
                        
                        overall_iter = overall_iter + 1;
                        seed_counter = seed_counter + 1;
                    end
                end
            end
        end
    end
    % End of city: record city elapsed time and update overall progress.
    city_elapsed = toc(city_tic);
    city_times(city_count) = city_elapsed;
    
    % Compute overall elapsed and estimated remaining time.
    elapsed_overall = toc(overall_tic);
    avg_iter_time = elapsed_overall / overall_iter;
    remaining_iter = total_iter - overall_iter;
    est_remaining = remaining_iter * avg_iter_time;
    
    % Print city summary on a new line.
    fprintf('City %s complete. Elapsed: %s | Estimated remaining: %s | Overall progress: %.1f%% (%d of %d iterations)\n', ...
        char(city_str), sec2time(city_elapsed), sec2time(est_remaining), (overall_iter/total_iter)*100, overall_iter, total_iter);
end

total_elapsed = toc(overall_tic);
fprintf('\nSimulation complete.\nTotal elapsed time: %s\n', sec2time(total_elapsed));
fprintf('Summary per city:\n');
for i = 1:numel(cities_str_array)
    fprintf('   %s: %s\n', char(cities_str_array(i)), sec2time(city_times(i)));
end

%% Helper function: sec2time converts seconds to HH:MM:SS format.
function time_str = sec2time(sec)
    hrs = floor(sec / 3600);
    mins = floor(mod(sec,3600) / 60);
    secs = mod(sec,60);
    time_str = sprintf('%02d:%02d:%05.2f', hrs, mins, secs);
end
