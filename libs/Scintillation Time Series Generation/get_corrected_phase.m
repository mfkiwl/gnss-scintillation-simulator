function [phase_int, n_interp] = get_corrected_phase(psi)
% get_corrected_phase - Construct an unwrapped phase from a complex field using FFT-based interpolation.
%
% Syntax:
%   [phase_int, n_interp, diff_w, wdiff] = get_corrected_phase(psi)
%
% Description:
%   This function computes the unwrapped phase of a complex field 'psi' by applying
%   FFT-based interpolation using MATLAB's built-in interpft function, along with the 
%   built-in angle and unwrap functions for phase processing. The interpolation factor is 
%   increased iteratively until the RMS difference between successive down-sampled phase 
%   values falls below 1 radian.
%
% Inputs:
%   psi - Complex field array, where each element is given as:
%         psi = |psi| * exp(1i * phase)
%
% Outputs:
%   phase_int - Unwrapped phase computed on the interpolated grid.
%
%   n_interp  - Interpolation factor used to achieve the error threshold.
%
%   diff_w    - Differences between successive wrapped phase samples.
%
%   wdiff     - Adjusted wrapped phase differences (mapped to [-pi, pi]).
%
% Example:
%   psi = exp(1i * linspace(0, 4*pi, 100));
%   [phase_int, n_interp, diff_w, wdiff] = get_corrected_phase(psi);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    max_interp = 10;     % Maximum allowed interpolation factor
    n_interp = 0;        
    error_val = Inf;     % Initialize RMS error as infinity
    prev_phase = [];
    
    nsamps = length(psi);
    
    while error_val > 1 && n_interp < max_interp
        n_interp = n_interp + 1;
        
        % Calculate the new length for interpolation.
        new_length = n_interp * nsamps;
        
        % FFT-based interpolation using MATLAB's built-in interpft.
        psi_interpolated = interpft(psi, new_length);

        % Compute the wrapped phase using the built-in angle function.
        wrapped_phase = angle(psi_interpolated);
        
        % Compute the unwrapped phase using the built-in unwrap function.
        phase_int_full = unwrap(wrapped_phase);
        
        % Down-sample to match the original sample positions.
        phase_original = phase_int_full(1:n_interp:end);
        
        % Compute RMS error if a previous phase exists.
        if ~isempty(prev_phase)
            error_val = rms(prev_phase - phase_original);
        end
        
        % Update the previous phase for the next iteration.
        prev_phase = phase_original;
    end
    
    % Final outputs:
    % phase_int is the unwrapped phase on the interpolated grid.
    phase_int = phase_int_full(1:n_interp:end);
end
