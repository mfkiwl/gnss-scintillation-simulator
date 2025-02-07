function norm_phase_sdf = get_norm_phase_sdf(mu, irr_param)
% get_norm_phase_sdf
%
% Syntax:
%   norm_phase_sdf = get_norm_phase_sdf(mu, irr_param)
%
% Description:
%   Compute the normalized phase spectral density function (SDF) using
%   a two-component power-law approach in a normalized wavenumber array 
%   defined by μ, based on parameters specified in `irr_param`. 
%   The function enforces a symmetric SDF around the 0 frequency by using 
%   abs(μ). The singularity at μ = 0 is explicitly set to zero to avoid 
%   undesirable distortions of the phase realizations.
%
% Inputs:
%   μ - Frequency array (e.g., spanning -50 to 50).
%
%   irr_param - Structure containing fields:
%       U   : Spectral strength
%       p1  : Exponent for the lower-frequency power-law region
%       p2  : Exponent for the higher-frequency power-law region
%       μ0  : Transition frequency between the two power-law regimes
%
% Outputs:
%   norm_phase_sdf - Vector of the same size as μ, containing the
%                    normalized phase SDF values.
%
% Notes:
%   - The absolute value of μ is used here to ensure that the SDF is
%     symmetric around the 0 frequency. Using the true (signed) μ as
%     proposed by equation (9) of [1] would result in an asymmetric SDF.
%
%   - The value at μ = 0 is set to zero (same approach used in the 
%     `GenPSRealization.m` script). This is feasible since the function 
%     is undefined at μ = 0 due to the power-law form, and setting it 
%     to zero prevents unwanted artifacts in the phase realizations.
%
% Example:
%   mu = -50:0.01:50;
%   irr_param.U   = 1;
%   irr_param.p1  = 2.0;
%   irr_param.p2  = 3.5;
%   irr_param.mu0 = 5;
%   sdf = get_norm_phase_sdf(mu, irr_param);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Preallocate the normalized phase spectral density function (SDF) values.
    norm_phase_sdf = zeros(size(mu));

    % Indices for below and above μ0
    idx_below = (abs(mu) <= irr_param.mu0);
    idx_above = (abs(mu) > irr_param.mu0);

    % First case
    norm_phase_sdf(idx_below) = ...
        irr_param.U .* abs(mu(idx_below)).^(-irr_param.p1);

    % Second case
    norm_phase_sdf(idx_above) = ...
        irr_param.U .* (irr_param.mu0^(irr_param.p2 - irr_param.p1)) ...
        .* abs(mu(idx_above)).^(-irr_param.p2);

    % Manually set the singular point at μ=0 to zero
    norm_phase_sdf(mu == 0) = 0;
end
