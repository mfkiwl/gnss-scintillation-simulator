function [extrapolated_irregularity_parameters, rhof_veff_vector] = ...
    freq_extrapolate_all(irregularity_params, general_parameters, rho_oveff_l1)
% freq_extrapolate_all Applies the frequency-extrapolation described by Eqs. (12), 
% (13), and (14) in [1] to 'Severe', 'Moderate', and 'Mild' irregularity parameters 
% for L2 and L5, using L1 as the reference.
%
% Syntax:
%   [extrapolated_irregularity_parameters, rhof_veff_vector] = ...
%       freq_extrapolate_all(irregularity_params, general_parameters, rho_oveff_l1)
%
% Description:
%   This function takes:
%       - irregularity_params, which must contain the fields:
%           .Severe, .Moderate, .Mild
%         each of these a struct with:
%             U   - Turbulence strength at L1,
%             mu0 - Normalized break wavenumber at L1,
%             p1  - Low-frequency spectral index,
%             p2  - High-frequency spectral index.
%       - general_parameters, a struct returned by get_general_parameters(),
%         providing the frequencies for [L1, L2, L5] in gps_bands.
%       - rho_oveff_l1, the ratio (rho_F / v_eff) known at L1.
%
%   The routine outputs:
%       1) extrapolated_irregularity_parameters, which stores newly computed
%          values of U and mu0 for each condition (Severe, Moderate, Mild)
%          at L2 and L5, while p1 and p2 remain unchanged.
%       2) rhof_veff_vector, a 1x3 vector containing (rho_F / v_eff) at L1, L2,
%          and L5, respectively, scaled according to Eq. (13) of [1].
%
%   For each condition, U and mu0 are extrapolated from L1 to L2 and L1 to L5
%   following the piecewise definitions from Eq. (14) and mu0 scaling from Eq. (12).
%
% Inputs:
%   irregularity_params     - Struct with fields .Severe, .Moderate, .Mild, each 
%                             containing U, mu0, p1, p2.
%   general_parameters      - Struct with at least .gps_bands, typically returned 
%                             by get_general_parameters().
%   rho_oveff_l1           - Ratio (rho_F / v_eff) at L1 (scalar).
%
% Outputs:
%   extrapolated_irregularity_parameters - Struct with two subfields: .L2 and .L5.
%       Each subfield has .Severe, .Moderate, and .Mild storing the updated 
%       U, mu0, p1, p2 for that frequency.
%
%   rhof_veff_vector       - 1x3 vector [rho_oveff_L1, rho_oveff_L2, rho_oveff_L5].
%
% Example:
%   irr_params = get_irregularity_parameters();
%   gen_params = get_general_parameters();
%   rho_oveff_l1 = 1.0;
%   [extrap_params, rho_vec] = freq_extrapolate_all(irr_params, gen_params, rho_oveff_l1);
%
% Notes: 
%   - This is an adaptation of the function `FreqExtrapolate.m` available
%   at the `gnss-scintillation-simulator-2-param` repository available at
%   https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param.
%   There were neither credits and a proper documentation on the original
%   code. However, as "Joy" commented on it, i'm assuming that Yu Jiao () have
%   participated on its development (https://www.researchgate.net/profile/Yu-Jiao-2).
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com
%
% [1] Xu D, Morton YTJ, Rino CL, Carrano CS, Jiao Y. 
%     A two-parameter multifrequency GPS signal simulator for strong equatorial 
%     ionospheric scintillation: modeling and parameter characterization. 
%     NAVIGATION. 2020; 67: 181–195. https://doi.org/10.1002/navi.350

    % Extract frequencies (L1, L2, L5) from general_parameters
    freq_l1 = general_parameters.gps_bands(1);
    freq_l2 = general_parameters.gps_bands(2);
    freq_l5 = general_parameters.gps_bands(3);
    
    extrapolated_irregularity_parameters = struct('L2', [], 'L5', []);
    
    % Scale the reference ratio (Eq. (13)) for L2 and L5
    rho_oveff_l2 = rho_oveff_l1 * sqrt(freq_l1 / freq_l2);
    rho_oveff_l5 = rho_oveff_l1 * sqrt(freq_l1 / freq_l5);
    rhof_veff_vector = [rho_oveff_l1, rho_oveff_l2, rho_oveff_l5];
    
    % Iterate over each condition (Severe, Moderate, Mild)
    for condition = {'Severe','Moderate','Mild'}
        cond_str = condition{1};
        U_f1   = irregularity_params.(cond_str).U;
        mu0_f1 = irregularity_params.(cond_str).mu0;
        p1     = irregularity_params.(cond_str).p1;
        p2     = irregularity_params.(cond_str).p2;
    
        % Extrapolate from L1 to L2 and from L1 to L5
        [U_L2, mu0_L2] = local_extrapolate(U_f1, mu0_f1, p1, p2, freq_l1, freq_l2);
        [U_L5, mu0_L5] = local_extrapolate(U_f1, mu0_f1, p1, p2, freq_l1, freq_l5);
    
        % Store results in the output struct
        extrapolated_irregularity_parameters.L2.(cond_str) = ...
            struct('U', U_L2, 'mu0', mu0_L2, 'p1', p1, 'p2', p2);
        extrapolated_irregularity_parameters.L5.(cond_str) = ...
            struct('U', U_L5, 'mu0', mu0_L5, 'p1', p1, 'p2', p2);
    end
end

function [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2)
    % It is important to comment here that some papers regarding the TPPSM
    % uses different notations on its definition.
    % 
    % Some papers defines the universal turbulance strength parameter `U` as:
    % U = Cpp * { 1,                  \mu_0 <= 1
    %           { \mu_0 ^{p_2 - p_1}, \mu_0 > 1
    % Which means that Cpp = U_1 defined in equation (5) of [1].
    % Therefore, we can calculate Cpp for the case when \mu_0 > 1 
    % by Cpp = U / (\mu_0 ^{p_2 - p_1}).
    %
    % The key idea here is that it is possible to calculate U_1=Cpp and U_2
    % by using a standard universal turbulance parameter `U` for simplicity.

    if mu0_f1 >= 1
        Cpp = U_f1; 
    else
        Cpp = U_f1 / mu0_f1^(p2 - p1); 
    end
    Cpp_f2 = Cpp * (f1 / f2)^(0.5 * p1 + 1.5);
    mu0_f2 = mu0_f1 * sqrt(f1 / f2);
    if mu0_f2 >= 1
        U_f2 = Cpp_f2;
    else
        U_f2 = Cpp_f2 * (mu0_f2)^(p2 - p1); 
    end
end
