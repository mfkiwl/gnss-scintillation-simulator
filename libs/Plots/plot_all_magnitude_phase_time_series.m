function plot_all_magnitude_phase_time_series(scint_field_struct, time_vector)
% plot_all_magnitude_phase
%
% Syntax:
%   plot_all_magnitude_phase(scint_field_struct, time_vector)
%
% Description:
%   Generates a **2x1 subplot** for each scintillation intensity level (Severe, Moderate, Weak).
%   - **First row**: Magnitude values (10*log10(abs(scint).^2)) for L1, L2, L5.
%   - **Second row**: Phase values (phase(scint)) for L1, L2, L5.
%   - **Titles include computed S4 values** for each frequency band.
%
% Inputs:
%   scint_field_struct - Struct containing scintillation realizations for each scenario:
%       .Severe.L1, .Severe.L2, .Severe.L5
%       .Moderate.L1, .Moderate.L2, .Moderate.L5
%       .Weak.L1, .Weak.L2, .Weak.L5
%
%   time_vector - 1D array containing the time values (same for all scenarios & frequencies).
%
% Outputs:
%   None (generates multiple figures and saves them as PDFs).
%
% Example:
%   plot_all_magnitude_phase(scint_field_struct, time_vector);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    scenarios = fieldnames(scint_field_struct);  % {'Severe', 'Moderate', 'Weak'}
    frequencies = {'L1', 'L2', 'L5'};
    colors = {'b', 'r', 'g'};  % L1 = Blue, L2 = Red, L5 = Green

    for i = 1:numel(scenarios)
        scenario = scenarios{i};
        
        % Create figure and set Helvetica as default font
        fig = figure('Name', sprintf('%s Magnitude & Phase', scenario), ...
                     'Position', [50, 50, 1400, 550], 'Color', 'w');
        set(fig, 'DefaultTextFontName', 'Helvetica');  
        tiledlayout(2, 1, "TileSpacing", "compact");

        % Initialize storage for legend labels and S4 values
        legend_labels = cell(1, numel(frequencies));
        s4_values = zeros(1, numel(frequencies));

        % Plot Magnitude (Row 1)
        nexttile;
        hold on;
        for j = 1:numel(frequencies)
            freq = frequencies{j};
            scint_field = scint_field_struct.(scenario).(freq);
            truncated_scint_field = scint_field(1:length(time_vector));
            % Compute S4
            s4_values(j) = get_S4(abs(truncated_scint_field).^2);
            
            plot(time_vector, 10*log10(abs(truncated_scint_field).^2), colors{j}, 'LineWidth', 1.2);
            legend_labels{j} = sprintf('%s (S_4=%.3f)', freq, s4_values(j));
        end
        xlabel('Time (s)', 'FontName', 'Helvetica');
        ylabel('Magnitude [dB]', 'FontName', 'Helvetica');
        title(sprintf('%s - Magnitude (Power in dB)', scenario), 'FontSize', 12, 'FontName', 'Helvetica');
        legend(legend_labels, 'Location', 'best', 'FontName', 'Helvetica');
        grid on;
        hold off;

        % Plot Phase (Row 2)
        nexttile;
        hold on;
        for j = 1:numel(frequencies)
            freq = frequencies{j};
            scint_field = scint_field_struct.(scenario).(freq);
            
            truncated_scint_field = scint_field(1:length(time_vector));
            plot(time_vector, phase(truncated_scint_field), colors{j}, 'LineWidth', 1.2);
        end
        xlabel('Time (s)', 'FontName', 'Helvetica');
        ylabel('Phase [rad]', 'FontName', 'Helvetica');
        title(sprintf('%s - Phase Evolution', scenario), 'FontSize', 12, 'FontName', 'Helvetica');
        legend(frequencies, 'Location', 'best', 'FontName', 'Helvetica');
        grid on;
        hold off;

        sgtitle(sprintf('Magnitude & Phase Analysis: %s', scenario), 'FontSize', 14, 'FontName', 'Helvetica');

        % Adjust figure paper size to match figure aspect ratio
        set(fig, 'PaperUnits', 'centimeters');
        fig_pos = get(fig, 'Position');
        fig_width_cm = fig_pos(3) * 2.54 / 96;  
        fig_height_cm = fig_pos(4) * 2.54 / 96;
        set(fig, 'PaperSize', [fig_width_cm fig_height_cm]);  
        set(fig, 'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

        % Save as PDF
        pdf_filename = sprintf('%s_Magnitude_Phase_Time_Series.pdf', scenario);
        print(fig, fullfile('..','figures',pdf_filename), '-dpdf', '-vector');
        
        disp(['Saved: ', fullfile('..','figures',pdf_filename)]);
    end
end
