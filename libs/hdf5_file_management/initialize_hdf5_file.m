function new_filename = initialize_hdf5_file(filename, city_names_str_array, attributes)
% initialize_hdf5_file
%
% Syntax:
%   new_filename = initialize_hdf5_file(filename, city_names_str_array, attributes)
%
% Description:
%   Checks if the HDF5 file already exists. If it does, the user is prompted to
%   either overwrite the entire file or create a new file version by appending a
%   suffix (e.g., _1, _2, ...). Then, a new HDF5 file is created under the
%   'hdf5_files' folder (three levels up) with one group for each city. In addition,
%   the provided attributes are attached to the root ("/") of the file.
%
% Inputs:
%   filename              - String. Name of the HDF5 file (e.g., 'simulation_results.h5').
%   city_names_str_array  - String array or cell array of city names.
%   attributes            - Struct containing attributes to attach to the root group.
%
% Outputs:
%   new_filename - String. The (possibly modified) filename used.
%
% Example:
%   attrs = struct('sim_date', '2014-01-02', 'simulation_time', 300);
%   new_filename = initialize_hdf5_file('simulation_results.h5', ["ascension_island", "hong_kong"], attrs);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Construct the full file path.
    file_path = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'hdf5_files', filename);

    % Check if the file exists.
    if exist(file_path, 'file')
        file_path_display = strrep(file_path, '\', '\\');
        prompt = sprintf('File %s exists. Overwrite the entire file (Y) or create a new version (N)? [Y/N]: ', file_path_display);
        user_choice = input(prompt, 's');
        if strcmpi(user_choice, 'Y')
            delete(file_path);
            new_filename = filename;
        else
            version = 1;
            [~, name, ext] = fileparts(filename);
            new_filename = sprintf('%s_%d%s', name, version, ext);
            new_file_path = fullfile(fileparts(file_path), new_filename);
            while exist(new_file_path, 'file')
                version = version + 1;
                new_filename = sprintf('%s_%d%s', name, version, ext);
                new_file_path = fullfile(fileparts(file_path), new_filename);
            end
            file_path = new_file_path;
        end
    else
        new_filename = filename;
    end

    % Create a new HDF5 file.
    file_id = H5F.create(file_path, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    
    % Create one group for each city.
    for i = 1:numel(city_names_str_array)
        % Convert city name to a character vector.
        city_name = char(city_names_str_array(i));
        group_path = ['/', city_name];
        group_id = H5G.create(file_id, group_path, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
        H5G.close(group_id);
    end
    
    H5F.close(file_id);
    
    % Attach attributes to the root group using high-level function.
    if nargin == 3 && ~isempty(attributes)
        attr_names = fieldnames(attributes);
        for i = 1:numel(attr_names)
            attr_name = attr_names{i};
            attr_value = attributes.(attr_name);
            h5writeatt(file_path, '/', attr_name, attr_value);
        end
    end

    fprintf('HDF5 file created: %s\n', file_path);
end
