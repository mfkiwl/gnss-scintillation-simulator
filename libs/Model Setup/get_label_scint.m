function label = get_label_scint(severity_str, scint_field, sampling_interval, varargin)
% get_label_scint
%
% Syntax:
%   label = get_label_scint(severity_str, scint_field, sampling_interval, varargin)
%
% Description:
%   This function computes the normalized autocorrelation of the intensity 
%   (I = |scint_field|^2) of a scintillation field at two specified lags defined 
%   by the 'lower_time' and 'upper_time' parameters. The autocorrelation values 
%   at these lags are compared against a threshold (default: exp(-1)) to classify 
%   the signal's decorrelation time as 'long', 'medium', or 'short'. The output 
%   label is constructed by concatenating the input severity category with the 
%   decorrelation classification (e.g., 'weak_long').
%
% Inputs:
%   severity_str      - String indicating the scintillation severity category.
%   scint_field       - Complex vector representing the scintillation field.
%   sampling_interval - Scalar specifying the time (in seconds) between samples.
%
% Optional Name-Value Pair Inputs:
%   'lower_time'      - Time (in seconds) corresponding to the lower lag for 
%                       autocorrelation evaluation (default: 0.5).
%   'upper_time'      - Time (in seconds) corresponding to the upper lag for 
%                       autocorrelation evaluation (default: 2.0).
%
% Outputs:
%   label - A string that combines the input severity category with the 
%           decorrelation classification ('_long', '_medium', or '_short').
%
% Error Conditions:
%   - Throws an error if lower_time is not less than upper_time.
%   - Throws an error if the signal is too short to compute autocorrelation at 
%     the specified lags.
%
% Notes:
%   - The autocorrelation is computed using MATLAB's autocorr function.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Parse optional parameters.
    p = inputParser;
    addParameter(p, 'lower_time', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'upper_time', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    lower_time = p.Results.lower_time;
    upper_time = p.Results.upper_time;
    
    % Ensure that lower_time is less than upper_time.
    if lower_time >= upper_time
        error('The lower_time must be less than the upper_time.');
    end
    
    % Compute the intensity of the scintillation signal.
    intensity = abs(scint_field).^2;
    N = length(intensity);
    
    % Calculate the lag indices corresponding to the lower and upper times.
    lag_lower = round(lower_time / sampling_interval);
    lag_upper = round(upper_time / sampling_interval);
    
    if lag_lower >= N || lag_upper >= N
        error('Signal is too short to compute autocorrelation at the specified lags.');
    end
    
    % Compute the autocorrelation using MATLAB's autocorr function.
    % We compute up to the upper lag.
    [acf, ~] = autocorr(double(intensity), 'NumLags', lag_upper);
    
    % Note: acf(1) corresponds to lag 0.
    R_lower = acf(lag_lower + 1);
    R_upper = acf(lag_upper + 1);
    
    % Define the threshold (default is exp(-1)).
    threshold = exp(-1);
    
    % Classify based on the threshold.
    if (R_lower > threshold) && (R_upper > threshold)
        decorrelation_char = 'long';
    elseif (R_lower > threshold) && (R_upper <= threshold)
        decorrelation_char = 'medium';
    else
        decorrelation_char = 'short';
    end
    
    % Combine severity string with the decorrelation classification.
    label = [char(severity_str), '_', decorrelation_char];
end
