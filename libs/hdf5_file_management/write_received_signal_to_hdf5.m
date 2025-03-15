function write_received_signal_to_hdf5(filename, city, data_real, data_imag, attributes, varargin)
% write_received_signal_to_hdf5
%
% Syntax:
%   write_received_signal_to_hdf5(filename, city, data_real, data_imag, attributes, varargin)
%
% Description:
%   Writes the received signal datasets (real and imaginary parts) to an existing
%   HDF5 file. The datasets are stored in the group corresponding to the receiver
%   location (city). Each dataset is a 30000x3 single array. All sweep parameters,
%   including PRN, drift velocity, CN0, severity, MC run, a [1x3] array of 
%   rhof_veff_ratios for L1, L2, and L5, and irregularity parameters (provided in 
%   the "irr_params" field of attributes) are attached to the dataset as attributes.
%   Optional parameters allow specifying the chunk size and deflate compression level.
%
% Inputs:
%   filename   - String. Name of the HDF5 file (e.g., 'simulation_results.h5').
%   city       - String. Name of the receiver city (must correspond to an existing group).
%   data_real  - 30000x3 single array containing the real part.
%   data_imag  - 30000x3 single array containing the imaginary part.
%   attributes - Struct. Must contain at least the following fields:
%                    prn              - Satellite PRN identifier.
%                    drift_velocity   - Numeric drift velocity (e.g., 25).
%                    cn0              - Numeric carrier-to-noise ratio in dB-Hz (e.g., 25).
%                    severity         - String. Severity level (e.g., 'weak').
%                    mc_run           - Numeric. Monte Carlo run number.
%                    rhof_veff_ratios - [1x3] double array with the rhof_veff ratios for L1, L2, and L5.
%                    irr_params       - Struct containing irregularity parameters (e.g., U, mu0, p1, p2).
%   varargin   - (Optional) Name-value pairs:
%                    'chunk_size'    - Numeric vector specifying the chunk size.
%                    'deflate_level' - Scalar numeric specifying the deflate compression level (0-9).
%
% Outputs:
%   (None) - The function writes datasets to the HDF5 file.
%
% Example:
%   attrs = struct('prn', 'prn_3', 'drift_velocity', 25, 'cn0', 25, 'severity', 'weak', 'mc_run', 1, ...
%                  'rhof_veff_ratios', [0.1 0.05 0.02], ...
%                  'irr_params', struct('U', 0.1, 'mu0', 0.05, 'p1', 0.001, 'p2', 0.002));
%   write_received_signal_to_hdf5('simulation_results.h5', 'ascension_island', data_real, data_imag, attrs, ...
%       'chunk_size', [10000, 3], 'deflate_level', 4);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Construct the full file path using the same folder structure.
    file_path = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'hdf5_files', filename);
    
    % Set up optional parameters.
    p = inputParser;
    addParameter(p, 'chunk_size', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'deflate_level', 0, @(x) isnumeric(x) && isscalar(x) && (x >= 0 && x <= 9));
    parse(p, varargin{:});
    chunk_size = p.Results.chunk_size;
    deflate_level = p.Results.deflate_level;
    
    % Build dataset base names using fields from attributes.
    base_name_real = sprintf('%s_vdrift_%d_cn0_%d_%s_mc%d_real', attributes.prn, attributes.drift_velocity, attributes.carrier_to_noise_ratio, attributes.severity, attributes.monte_carlo_run);
    base_name_imag = sprintf('%s_vdrift_%d_cn0_%d_%s_mc%d_imag', attributes.prn, attributes.drift_velocity, attributes.carrier_to_noise_ratio, attributes.severity, attributes.monte_carlo_run);
    
    % Construct full dataset paths (city group should already exist).
    full_path_real = sprintf('/%s/%s', char(city), base_name_real);
    full_path_imag = sprintf('/%s/%s', char(city), base_name_imag);
    
    % Create the datasets with optional chunking and compression.
    if ~isempty(chunk_size)
        h5create(file_path, full_path_real, size(data_real), 'Datatype', class(data_real), 'ChunkSize', chunk_size, 'Deflate', deflate_level);
        h5create(file_path, full_path_imag, size(data_imag), 'Datatype', class(data_imag), 'ChunkSize', chunk_size, 'Deflate', deflate_level);
    else
        h5create(file_path, full_path_real, size(data_real), 'Datatype', class(data_real), 'Deflate', deflate_level);
        h5create(file_path, full_path_imag, size(data_imag), 'Datatype', class(data_imag), 'Deflate', deflate_level);
    end
    
    % Write the data to the datasets.
    h5write(file_path, full_path_real, data_real);
    h5write(file_path, full_path_imag, data_imag);
    
    % Loop over all fields in the attributes struct and attach them to both datasets.
    att_fields = fieldnames(attributes);
    for i = 1:length(att_fields)
        att_name = att_fields{i};
        att_value = attributes.(att_name);
        h5writeatt(file_path, full_path_real, att_name, att_value);
        h5writeatt(file_path, full_path_imag, att_name, att_value);
    end
    
    %fprintf('Datasets written: %s and %s in group /%s in file %s\n', full_path_real, full_path_imag, char(city), file_path);
end
