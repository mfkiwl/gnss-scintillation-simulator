function detrended_phase_realization = get_phase_realization(norm_phase_sdf, D_mu, nfft, seed)
% get_phase_realization
%
% Syntax:
%   phase_realization = get_phase_realization(norm_phase_sdf, D_mu, seed)
%
% Description:
%   This function returns the real part of the inverse Fourier transform of the 
%   product between a square root of the normalized phase spectral density function 
%   (norm_phase_sdf) and a Gaussian random complex vector (xi). It uses a fixed seed 
%   for reproducibility, ensuring the same results each time the function is called.
%
% Inputs:
%   norm_phase_sdf - Normalized phase spectral density function values
%   D_mu           - Vector of frequency spacing or related DFT increment
%   seed           - Seed for the random number generator to ensure repeatability
%
% Outputs:
%   phase_realization - The real-valued phase realization in the time domain
%
% Notes:
%   - The variance of xi is sqrt(2). In the literature, it is recommended that xi 
%     should have variance 1, be Hermitian, and be white. However, the following 
%     codes:
%
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/master/Libraries/GenScintFieldRealization/GenPSRealization.m
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator/blob/master/matlab/PropCodes/turbsim1.m
%
%     do not enforce variance = 1 or the Hermitian property. On the
%     contrary, they introduces a scaling factor of 1.22 (on the first code) 
%     and \sqrt(2) (on the second code), respectively on each of the phase 
%     realizations, which do not make any physical sense, without 
%     commenting over this approach. It is important to comment also that
%     these scaling factors acts altering the value of the turbulance
%     strength Cpp significantly, which impacts the whole reliability of
%     the code.
%
%   - Using the Hermitian property as indicated in the TPPSM literature caused 
%     distortions at the beginning and end of the final scintillation time series. 
%     Hence, it might not be strictly necessary to use the Hermitian property.
%
%   - The expression `fftshift(fft(fftshift(root_norm_phase_sdf .* xi)))` produces 
%     a complex time series with nonzero real and imaginary parts if xi even if is 
%     Hermitian. This conflicts with literature that suggests the phase realization 
%     would be purely real if xi were Hermitian [2, Equations 2.39 - 2.43].
%
%   - The only way found to match the periodogram and the value of S4 to the 
%     theoretical spectrum provided by Ispectrum.exe was to use xi with variance 
%     sqrt(2). This may be due to "power leakage" in 
%     `fftshift(fft(fftshift(root_norm_phase_sdf .* xi)))`, 
%     where the imaginary part can be significant.
%
% Example:
%   % Example usage:
%   norm_phase_sdf = rand(1, 256);  % Example vector
%   D_mu = 0.1;                     % Example spacing
%   seed = 12345;                   % Example seed
%   phase_realization = get_phase_realization(norm_phase_sdf, D_mu, seed);
%
% References:
% [1] Rino C, Breitsch B, Morton Y, Jiao Y, Xu D, Carrano C. A compact 
%     multi-frequency GNSS scintillation model. NAVIGATION. 2018; 65: 
%     563–569. https://doi.org/10.1002/navi.263
% [2] C. Rino, The Theory of Scintillation with Applications in Remote
%     Sensing. John Wiley & Sons, 2011
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Generate Gaussian random complex vector (xi). 
    % Note that the variance of xi is sqrt(2).
    xi = randn(1, nfft) + 1i * randn(1, nfft);

    % Compute the square-root of the normalized phase SDF.
    root_norm_phase_sdf = sqrt(norm_phase_sdf * D_mu / (2*pi));

    % Obtain the phase realization by taking the real part of the
    % inverse Fourier transform (with shift operations), which corresponds to equation (22) of (1).
    phase_realization = real(fftshift(fft(fftshift(root_norm_phase_sdf .* xi))));

    %Remove linear trend to force segment too segment continuity
    % NOTE (Rodrigo): I've let this part the same as it is done in the
    % TPPSM original code:
    % https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/master/Libraries/GenScintFieldRealization/GenScintFieldRealization.m
    linear_trend = linex( ...
        1:nfft, ...
        1, ...
        nfft, ...
        phase_realization(1), ...
        phase_realization(nfft));
    detrended_phase_realization = phase_realization - linear_trend;
end
