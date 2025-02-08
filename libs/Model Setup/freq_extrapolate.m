function [extrapolated_irr_params, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params, gen_params, rho_oveff_ratio_l1)
% freq_extrapolate
%
% Syntax:
%   [extrapolated_irr_params, rhof_veff_ratio_vector] = ...
%       freq_extrapolate(irr_params, gen_params, rho_oveff_ratio_l1)
%
% Description:
%   This function extrapolates the irregularity parameters (U, mu0, p1, p2) 
%   from a reference frequency (L1) to two other frequencies (L2 and L5). 
%   The scaling follows the conventions in [1]:
%       - Equation (12) for scaling the break wavenumber mu0.
%       - Equation (13) for scaling the ratio (rho_F / v_eff).
%       - Equation (14) for the piecewise definition of the universal 
%         turbulence strength parameter U.
%
%   Specifically:
%   1) mu0 at L2 and L5 is obtained by multiplying the L1 value by sqrt(L1_frequency / target_frequency) [Eq. (12)].
%   2) (rho_F / v_eff) for L2 and L5 is determined by multiplying the L1 value 
%      by sqrt(L1_frequency / target_frequency) [Eq. (13)].
%   3) U is computed via the four-case piecewise expression [Eq. (14)], depending 
%      on whether mu0 >= 1 at the reference or target frequency.
%
% Inputs:
%   irr_params - Struct with fields:
%       .U   - Turbulence strength at L1
%       .mu0 - Normalized break wavenumber at L1
%       .p1  - Low-frequency spectral index
%       .p2  - High-frequency spectral index
%
%   gen_params  - Struct returned by get_gen_params(), containing 
%                        the field .gps_bands, a 1x3 vector with the frequencies 
%                        [L1, L2, L5].
%
%   rho_oveff_ratio_l1  - Ratio (rho_F / v_eff) known at L1 (scalar).
%
% Outputs:
%   extrapolated_irr_params - Struct with fields .L2 and .L5,
%       each containing .U, .mu0, .p1, and .p2. These represent the 
%       extrapolated parameters at the corresponding frequency.
%
%   rhof_veff_ratio_vector - 1x3 vector [rho_oveff_L1, rho_oveff_L2, rho_oveff_L5].
%
% Notes:
%   - This code is adapted from the function FreqExtrapolate.m available at:
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/master/...
%       Libraries/GenScintFieldRealization/FreqExtrapolate.m
%
% Example:
%   % Suppose we have a single set of irregularity parameters:
%   irr_params.U   = 1.5;
%   irr_params.mu0 = 0.8;
%   irr_params.p1  = 2.0;
%   irr_params.p2  = 3.5;
%
%   % Get the general parameters (which includes gps_bands for L1, L2, L5)
%   gen_params = get_gen_params();
%
%   % Known ratio (rho_F / v_eff) at L1
%   rho_oveff_l1 = 1.0;
%
%   % Extrapolate
%   [extrap_params, rho_vec] = freq_extrapolate(irr_params, gen_params, rho_oveff_l1);
%
% References:
%   [1] Xu D, Morton YTJ, Rino CL, Carrano CS, Jiao Y. 
%       "A two-parameter multifrequency GPS signal simulator for strong 
%        equatorial ionospheric scintillation: modeling and parameter 
%        characterization." NAVIGATION. 2020; 67: 181–195. 
%       https://doi.org/10.1002/navi.350
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Extract frequencies (L1, L2, L5) from gen_params.
    freq_l1 = gen_params.gps_bands(1);
    freq_l2 = gen_params.gps_bands(2);
    freq_l5 = gen_params.gps_bands(3);

    % Scale the reference ratio (rho_F / v_eff) for L2 and L5 [Eq. (13)].
    rho_oveff_ratio_l2 = rho_oveff_ratio_l1 * sqrt(freq_l1 / freq_l2);
    rho_oveff_ratio_l5 = rho_oveff_ratio_l1 * sqrt(freq_l1 / freq_l5);
    rhof_veff_ratio_vector = [rho_oveff_ratio_l1, rho_oveff_ratio_l2, rho_oveff_ratio_l5];

    % Extract parameters at L1.
    U_f1   = irr_params.U;
    mu0_f1 = irr_params.mu0;
    p1     = irr_params.p1;
    p2     = irr_params.p2;

    % Extrapolate from L1 to L2 and from L1 to L5.
    [U_L2, mu0_L2] = local_extrapolate(U_f1, mu0_f1, p1, p2, freq_l1, freq_l2);
    [U_L5, mu0_L5] = local_extrapolate(U_f1, mu0_f1, p1, p2, freq_l1, freq_l5);

    % Store results in the output struct.
    extrapolated_irr_params = struct('L2', [], 'L5', []);
    extrapolated_irr_params.L2 = struct('U', U_L2,  ...
                                                     'mu0', mu0_L2, ...
                                                     'p1', p1, ...
                                                     'p2', p2);
    extrapolated_irr_params.L5 = struct('U', U_L5,  ...
                                                     'mu0', mu0_L5, ...
                                                     'p1', p1, ...
                                                     'p2', p2);

end

function [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2)
% local_extrapolate
%
% Syntax:
%   [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2)
%
% Description:
%   Applies Eq. (12) from [1] to scale mu0 between frequencies, and applies 
%   the piecewise expression of Eq. (14) to determine the new turbulence 
%   strength U. The specific case depends on whether mu0 >= 1 at the 
%   reference frequency (f1) or at the target frequency (f2).
%
%   - mu0 is scaled by: mu0_f2 = mu0_f1 * sqrt(f1 / f2)    [Eq. (12)]
%   - U is scaled by a factor of (f1/f2)^exponent, where 
%     exponent depends on p1 or p2, depending on the four cases in [Eq. (14)].
%
% Inputs:
%   U_f1   - Turbulence strength at the reference frequency f1
%   mu0_f1 - Break wavenumber at the reference frequency f1
%   p1     - Low-frequency spectral index
%   p2     - High-frequency spectral index
%   f1     - Reference frequency
%   f2     - Target frequency
%
% Outputs:
%   U_f2   - Extrapolated turbulence strength at the target frequency f2
%   mu0_f2 - Extrapolated break wavenumber at the target frequency f2
%
% Example:
%   U_f1   = 1.5;
%   mu0_f1 = 0.8;
%   p1     = 2.0;
%   p2     = 3.5;
%   f1     = 1.57542e9;  % Example L1 frequency
%   f2     = 1.22760e9;  % Example L2 frequency
%   [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Exponents for (f1/f2) based on p1 or p2 [Eq. (14)].
    exponent_p1 = 0.5 * p1 + 1.5;
    exponent_p2 = 0.5 * p2 + 1.5;

    % Scale mu0 from f1 to f2 [Eq. (12)].
    mu0_f2 = mu0_f1 * sqrt(f1 / f2);

    % Determine the cases for Eq. (14).
    cond_f1 = (mu0_f1 >= 1);
    cond_f2 = (mu0_f2 >= 1);

    if  cond_f1 &&  cond_f2
        % Case 1: mu0(f1) >= 1, mu0(f2) >= 1
        U_f2 = U_f1 * (f1/f2)^exponent_p1;
    elseif ~cond_f1 &&  cond_f2
        % Case 2: mu0(f1) < 1,  mu0(f2) >= 1
        U_f2 = U_f1 * (mu0_f1)^(p1 - p2) * (f1/f2)^exponent_p1;
    elseif  cond_f1 && ~cond_f2
        % Case 3: mu0(f1) >= 1, mu0(f2) < 1
        U_f2 = U_f1 * (mu0_f2)^(p2 - p1) * (f1/f2)^exponent_p1;
    else
        % Case 4: mu0(f1) < 1,  mu0(f2) < 1
        U_f2 = U_f1 * (f1/f2)^exponent_p2;
    end
end
