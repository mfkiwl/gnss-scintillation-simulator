function label = get_label_scint(scint_field, sampling_interval, varargin)
% get_label_scint
%
% Syntax:
%   label = get_label_scint(scint_field, sampling_interval, varargin)
%
% Description:
%   Computes the normalized autocorrelation of the signal intensity of the scintillation
%   field at two specified time lags ("lower" and "upper") and assigns a label based on these
%   values.
%
%   The intensity of the scintillation signal is defined as:
%       I(n) = |scint_field(n)|^2
%
%   The autocorrelation function (ACF) of I(n) at lag m is defined as:
%       R(m) = (sum_{n=1}^{N-m} I(n) * I(n+m)) / (sum_{n=1}^{N} I(n)^2)
%
%   This normalization ensures that R(0) = 1, which is equivalent to the MATLAB command:
%       xcorr(I, 'coeff')
%
%   Our approach computes R(m) directly at the two desired lags (lower and upper). The 
%   function then assigns a label according to the following criteria:
%       - 'long'   if R(lower) > threshold and R(upper) > threshold,
%       - 'medium' if R(lower) > threshold but R(upper) < threshold,
%       - 'short'  otherwise.
%
% Inputs:
%   scint_field       - Complex vector representing the scintillation field.
%   sampling_interval - Scalar, time (in seconds) between samples.
%
% Optional Name-Value Pair Inputs:
%   'lower_time'      - Time (in seconds) for the lower lag evaluation (default: 0.5).
%   'upper_time'      - Time (in seconds) for the upper lag evaluation (default: 2.0).
%
% Outputs:
%   label - String: 'long', 'medium', or 'short'.
%
% Example:
%   field = randn(30000,1) + 1j*randn(30000,1);
%   label = get_label_scint(field, 0.01, 'lower_time', 0.5, 'upper_time', 2.0);
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
    
    % Compute the intensity of the scintillation signal.
    intensity = abs(scint_field).^2;
    N = length(intensity);
    
    % Calculate the lag indices corresponding to the lower and upper times.
    lag_lower = round(lower_time / sampling_interval);
    lag_upper = round(upper_time / sampling_interval);
    
    if lag_lower >= N || lag_upper >= N
        error('Signal is too short to compute autocorrelation at the specified lags.');
    end
    
    % Compute the normalized autocorrelation at the specified lags:
    % R(m) = (sum_{n=1}^{N-m} I(n)*I(n+m)) / (sum_{n=1}^{N} I(n)^2)
    R_lower = sum(intensity(1:N-lag_lower) .* intensity(1+lag_lower:N)) / sum(intensity.^2);
    R_upper = sum(intensity(1:N-lag_upper) .* intensity(1+lag_upper:N)) / sum(intensity.^2);
    
    % Define a threshold (default is exp(-1)).
    threshold = exp(-1);
    
    % Assign label based on the computed autocorrelation values:
    if (R_lower > threshold) && (R_upper > threshold)
        label = 'long';
    elseif (R_lower > threshold) && (R_upper < threshold)
        label = 'medium';
    else
        label = 'short';
    end
end
