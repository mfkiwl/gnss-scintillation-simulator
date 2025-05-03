function phase_int = get_corrected_phase(psi)
% get_corrected_phase 
% 
% Construct an unwrapped phase from a complex field using FFT-based interpolation.
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
%   psi - Complex field that represents the scintillation effects observed
%         at the receiver's plane.
%
% Outputs:
%   phase_int - Unwrapped phase computed on the interpolated grid.
%
% Example:
%   phase_int = get_corrected_phase(psi);
%
% Notes:
%   - This function was deeply inspired on the `ConstructPhase.m` provided
%   by the original code developed by the Colardo Boulder University group,
%   which can be accessed in the following links:
%     - https://github.com/ita-gnss-lab/gnss-scintillation-simulator/blob/ecbff9cd1b8bb32101a43e447fede7196128c72b/matlab/PropCodes/ConstructPhase.m
%     - https://github.com/ita-gnss-lab/gnss-scintillation-simulator/blob/ecbff9cd1b8bb32101a43e447fede7196128c72b/matlab/PropCodes/fftInterp.m
%     - https://github.com/ita-gnss-lab/gnss-scintillation-simulator/blob/ecbff9cd1b8bb32101a43e447fede7196128c72b/matlab/PropCodes/wrapPhase.m
%     - https://github.com/ita-gnss-lab/gnss-scintillation-simulator/blob/ecbff9cd1b8bb32101a43e447fede7196128c72b/matlab/PropCodes/unwrapPhase.m
%     - https://github.com/ita-gnss-lab/gnss-scintillation-simulator/blob/ecbff9cd1b8bb32101a43e447fede7196128c72b/matlab/Utilities/nicefftnum.m
%   - The functions `fftInterp`, `wrapPhase`, `unwrapPhase` and
%     `nicefftnum` were substituted by bult-in matlab functions, aiming to
%     enhance the manteinability of the code.
%   - For more details on the reasoning behind the approach used herein for
%   phase unwrapping correction, please refer to [1].
%
% References:
%   [1] Rino C, Breitsch B, Morton Y, Xu D, Carrano C. GNSS signal phase, 
%       TEC, and phase unwrapping errors. NAVIGATION. 2020; 67: 865–873. 
%       https://doi.org/10.1002/navi.396
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Maximum allowed interpolation factor
    % NOTE: The value was configured arbitrarily as 10, given that most of
    % the uncorrected phase unwraps for severe scintillation intensity can
    % be corrected with few iterations for strong scintillation events.
    max_interp = 10;
    % Initial interpolation step
    n_interp = 0;
    % Initialize RMS error as infinity
    error_val = Inf;
    % Initialize the previous phase time series vector
    prev_phase = [];
    % Amount of samples on the scintillation complex field.
    nsamps = length(psi);

    while error_val > 1 && n_interp < max_interp
        % Increment the ratio of interpolated samples.
        n_interp = n_interp + 1;
        
        % Calculate the new length for interpolation.
        new_length = n_interp * nsamps;
        
        % FFT-based interpolation using MATLAB's built-in `interpft` function.
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
    
    % phase_int is the unwrapped phase on the interpolated grid.
    phase_int = phase_int_full(1:n_interp:end);
end
