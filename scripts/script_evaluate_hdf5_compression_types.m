% script_evaluate_hdf5_compression_types.m
%
% Description:
% This script generates a scintillation time series for L1 (Severe case),
% converts it to single precision, and then saves it to an HDF5 file using
% various compression options. Since HDF5 does not natively support complex
% numbers, the real and imaginary parts are stored separately.
%
% Script Sections:
%
% 1. Add Required Paths
%    - Adds the necessary library and cache directories to the MATLAB path.
%
% 2. Model Setup & Time Series Generation (Severe, L1)
%    - Retrieves general simulation parameters and the rhof/veff ratio for L1.
%    - Retrieves the irregularity parameters and extrapolates them for the Severe case.
%    - Generates a complex double scintillation time series for L1 using the custom algorithm.
%
% 3. Convert to Single Precision
%    - Converts the generated complex double time series to single precision 
%      to facilitate HDF5 storage.
%
% 4. HDF5 Compression Testing
%    - Defines a function (test_hdf5_compression_complex) that:
%         * Splits the complex data into real and imaginary parts.
%         * Saves these parts to an HDF5 file using specified compression options.
%         * Measures and returns the save time, load time, and file size.
%
% 5. Define and Evaluate Compression Options
%    - Specifies various HDF5 compression options (with and without the Shuffle filter).
%    - Evaluates each option by saving and loading the data, and collects performance metrics.
%
% 6. Display Results
%    - Outputs the compression performance metrics (save time, load time, file size)
%      for each compression option.
%    - Identifies and displays the best options based on file size and load speed.
%
% Dependencies:
% This script relies on the following custom functions from the developed library:
%   - get_general_parameters
%   - get_rhof_veff_ratio
%   - get_irregularity_parameters
%   - freq_extrapolate
%   - get_scintillation_time_series
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

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

%% Convert to Single Precision (for HDF5 storage)
scint_L1_single = single(scint_L1);

%% Function to Test HDF5 Compression Performance for Complex Data
function [save_time, load_time, file_size] = test_hdf5_compression_complex(data, filename, deflateLevel, useShuffle)
    % data: complex single data
    % filename: output filename for HDF5 file
    % deflateLevel: integer 0-9 (0 means no compression)
    % useShuffle: boolean flag for using the Shuffle filter
    
    % Delete file if it exists, to avoid group/dataset conflicts.
    if exist(filename, 'file')
        delete(filename);
    end

    % Split data into real and imaginary parts
    real_data = real(data);
    imag_data = imag(data);
    dataSize = size(real_data);
    
    % Choose a chunk size that does not exceed the dataset dimensions.
    chunkRows = min(1000, dataSize(1));
    chunkSize = [chunkRows, dataSize(2)];
    
    % Build options list for h5create.
    opts = {};
    if useShuffle
        opts = [opts, {'Shuffle', true}];
    end
    if deflateLevel > 0
        opts = [opts, {'Deflate', deflateLevel}];
    end
    opts = [opts, {'ChunkSize', chunkSize}];
    
    % Save real and imaginary parts separately and measure time.
    tic;
    % Save real part
    h5create(filename, '/data/real', dataSize, 'Datatype', 'single', opts{:});
    h5write(filename, '/data/real', real_data);
    % Save imaginary part
    h5create(filename, '/data/imag', dataSize, 'Datatype', 'single', opts{:});
    h5write(filename, '/data/imag', imag_data);
    save_time = toc;
    
    % Measure load time: read both parts and reconstruct complex signal.
    tic;
    loaded_real = h5read(filename, '/data/real');
    loaded_imag = h5read(filename, '/data/imag');
    loaded_data = loaded_real + 1i * loaded_imag; %#ok<NASGU>
    load_time = toc;
    
    % Get file size in MB (the total file size)
    file_info = dir(filename);
    file_size = file_info.bytes / (1024^2);
end

%% Define Compression Options to Test
% Each option is defined by a label, a deflate level, and a flag for using Shuffle.
compressionOptions = {...
    struct('label', 'No_Compression', 'deflate', 0, 'shuffle', false), ...
    struct('label', 'Gzip_Level_1',    'deflate', 1, 'shuffle', false), ...
    struct('label', 'Gzip_Level_5',    'deflate', 5, 'shuffle', false), ...
    struct('label', 'Gzip_Level_9',    'deflate', 9, 'shuffle', false), ...
    struct('label', 'Gzip_Level_9_Shuffle', 'deflate', 9, 'shuffle', true) ...
};

%% Evaluate Compression Options
results = struct();
for k = 1:length(compressionOptions)
    opt = compressionOptions{k};
    filename = sprintf('test_single_%s.h5', strrep(opt.label, ' ', '_'));
    [sTime, lTime, fSize] = test_hdf5_compression_complex(scint_L1_single, filename, opt.deflate, opt.shuffle);
    results.(strrep(opt.label, ' ', '_')) = struct(...
        'save_time', sTime, ...
        'load_time', lTime, ...
        'file_size_MB', fSize ...
    );
end

%% Display Results
disp('=== HDF5 Compression Performance for Single Precision Complex Data ===');
optionLabels = fieldnames(results);
for i = 1:numel(optionLabels)
    optRes = results.(optionLabels{i});
    fprintf('%s:\n  Save Time: %.4f s\n  Load Time: %.4f s\n  File Size: %.2f MB\n\n', ...
        optionLabels{i}, optRes.save_time, optRes.load_time, optRes.file_size_MB);
end

fprintf('\nBest option based on file size:\n');
minSize = inf;
bestOption = '';
for i = 1:numel(optionLabels)
    sz = results.(optionLabels{i}).file_size_MB;
    if sz < minSize
        minSize = sz;
        bestOption = optionLabels{i};
    end
end
fprintf('  %s with %.2f MB\n', bestOption, minSize);

fprintf('\nBest option based on load speed:\n');
minLoad = inf;
bestOption = '';
for i = 1:numel(optionLabels)
    lt = results.(optionLabels{i}).load_time;
    if lt < minLoad
        minLoad = lt;
        bestOption = optionLabels{i};
    end
end
fprintf('  %s with %.4f s\n', bestOption, minLoad);
