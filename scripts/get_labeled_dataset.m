%% Initialization
clearvars; clc;

addpath(genpath(fullfile('..','libs')));
addpath(fullfile('..','cache'));

%% Parameter Setup
% Amount of drift velocities linearly sampled from 25 to 125m/s
drift_velocities_amount = 1;
% Amount of carrier-to-noise ratio values varying from 25 to 45 dB-Hz
carrier_to_noise_ratios_amount = 1;

monte_carlo_runs = 1;
% Mask angle for the line-of-sight satellites
mask_angle_deg = 15;

data_set_params = get_dataset_params(drift_velocities_amount, carrier_to_noise_ratios_amount, mask_angle_deg);

cities_str_array = string(fieldnames(data_set_params.specific));

severity_str_array = string(fieldnames(data_set_params.general.irr_params));
receiver_sampling_frequency = 2e7;
B = receiver_sampling_frequency/2;
rx_mean_power = 1;

frequency_str_array = ["L1", "L2", "L5"];

%% Main loop

for city_str = cities_str_array.'
    for prn_str = data_set_params.specific.(char(city_str)).prn_string_array
        for drift_idx = 1:drift_velocities_amount
            for severity_str = severity_str_array.'
                for carrier_to_noise_ratio = data_set_params.general.carrier_to_noise_ratios
                    for mc_run = 1 : monte_carlo_runs
                        real_received_signal_single_array = zeros(data_set_params.general.simulation_time / data_set_params.general.dt,1);
                        imag_received_signal_single_array = zeros(data_set_params.general.simulation_time / data_set_params.general.dt,1);
                        
                        for frequency_str = frequency_str_array

                            frequency_idx = find(frequency_str_array == frequency_str);

                            thermal_noise_single_array = get_thermal_noise( ...
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

                            scint_single_array = get_scintillation_time_series( ...
                                data_set_params.general, ...
                                irr_params, ...
                                rhof_veff_ratio, ...
                                'data_type', 'single');
                            
                            scint_single_array_sliced = scint_single_array(1:data_set_params.general.simulation_time/data_set_params.general.dt).';

                            received_signal_single_array = scint_single_array_sliced + thermal_noise_single_array;

                            real_received_signal_single_array(:,frequency_idx) = real(received_signal_single_array);
                            imag_received_signal_single_array(:,frequency_idx) = imag(received_signal_single_array);
                        end
                    end
                end
            end
        end
    end
end