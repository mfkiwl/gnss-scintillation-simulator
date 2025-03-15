function irregularity_params = get_irregularity_parameters(varargin)
% get_irregularity_parameters Returns a struct with irregularity parameters
% for different turbulence conditions.
%
% Syntax:
%   irregularity_params = get_irregularity_parameters()
%   irregularity_params = get_irregularity_parameters('Strong', severeParams, ...
%                                             'Moderate', moderateParams, ...
%                                             'Weak', WeakParams)
%
% Description:
%   This function outputs a struct (irregularity_params) with three fields: 'Strong',
%   'Moderate', and 'Weak'. Each field is itself a struct containing the
%   following irregularity parameters:
%       U   - Universal Turbulence Strength
%       mu0 - Normalized Break Wavenumber
%       p1  - Low-frequency Spectral Index
%       p2  - High-frequency Spectral Index
%
%   If no input is provided, default values are used for each condition.
%   Alternatively, the user can supply custom parameter values via name-value
%   pairs where the value is a struct with fields 'U', 'mu0', 'p1', and 'p2'.
%
% Inputs:
%   'Strong'   - (Optional) A struct with the fields:
%                    U   - Universal Turbulence Strength
%                    mu0 - Normalized Break Wavenumber
%                    p1  - Low-frequency Spectral Index
%                    p2  - High-frequency Spectral Index
%                Default: struct('U', 1.0, 'mu0', 1.0, 'p1', 3.0, 'p2', 4.0)
%
%   'Moderate' - (Optional) A struct with the same fields as above.
%                Default: struct('U', 0.8, 'mu0', 0.8, 'p1', 2.8, 'p2', 3.8)
%
%   'Weak'     - (Optional) A struct with the same fields as above.
%                Default: struct('U', 0.6, 'mu0', 0.6, 'p1', 2.5, 'p2', 3.5)
%
% Outputs:
%   irregularity_params - A struct with three fields ('Strong', 'Moderate', 'Weak'), each
%                containing a struct with irregularity parameters.
%
% Example:
%   % Use default parameters:
%   irregularity_params = get_irregularity_parameters();
%
%   % Override the Strong condition parameters:
%   customSevere = struct('U', 1.2, 'mu0', 1.1, 'p1', 3.2, 'p2', 4.1);
%   irregularity_params = get_irregularity_parameters('Strong', customSevere);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    %% Define default parameter structures for each turbulence condition
    defaultSevere   = struct('U', 2, 'mu0', 0.55, 'p1', 2.45, 'p2', 3.7);
    % TODO: Confirm wether we should adopt the same mu0, p1 and p2 for both
    % severe and moderate scintillation cases.
    defaultModerate = struct('U', 0.4, 'mu0', 0.7, 'p1', 2.7, 'p2', 3.3);
    % Note: The `defaultWeak` parameters is a one-component power law, since
    % p1 = p2. It is noteworthy to comment that mu0 have no impact on the
    % model in this case.
    defaultWeak     = struct('U', 0.05, 'mu0', 1, 'p1', 3, 'p2', 3);

    %% Create input parser and add optional parameters
    pInput = inputParser;
    
    % Validation function to check that the input is a 
    % struct with the required fields.
    validateConditionStruct = @(x) isstruct(x) && ...
        all(isfield(x, {'U', 'mu0', 'p1', 'p2'}));
    
    % Add parameters for each condition using the default structs
    addParameter(pInput, 'strong', defaultSevere, validateConditionStruct);
    addParameter(pInput, 'moderate', defaultModerate, validateConditionStruct);
    addParameter(pInput, 'weak', defaultWeak, validateConditionStruct);
    
    % Parse user-supplied name-value pairs
    parse(pInput, varargin{:});
    
    %% Build the output struct with the irregularity parameters for each condition
    irregularity_params = struct();
    irregularity_params.strong   = pInput.Results.strong;
    irregularity_params.moderate = pInput.Results.moderate;
    irregularity_params.weak     = pInput.Results.weak;
end
