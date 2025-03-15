function set_city_group_attributes(filename, city, attributes)
% set_city_group_attributes
%
% Syntax:
%   set_city_group_attributes(filename, city, attributes)
%
% Description:
%   Attaches the provided attributes to the group corresponding to a given city in the HDF5 file.
%
% Inputs:
%   filename   - String. Name of the HDF5 file (e.g., 'simulation_results.h5').
%   city       - String. Name of the city (must correspond to an existing group).
%   attributes - Struct. Fields represent attribute names and values.
%
% Example:
%   attrs = struct('location', 'tropical', 'altitude', 728);
%   set_city_group_attributes('simulation_results.h5', 'ascension_island', attrs);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Construct full file path using the same logic.
    file_path = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'hdf5_files', filename);
    
    % Construct the group path (ensure city is a char vector).
    group_path = sprintf('/%s', char(city));
    
    % Loop over attribute fields and attach them.
    attr_names = fieldnames(attributes);
    for i = 1:numel(attr_names)
        attr_name = attr_names{i};
        attr_value = attributes.(attr_name);
        h5writeatt(file_path, group_path, attr_name, attr_value);
    end
    
    %fprintf('Attributes set for group %s in file %s\n', group_path, file_path);
end
