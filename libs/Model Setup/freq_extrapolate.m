function [extrapolated_irr_params, rhof_veff_ratio_vector] = ...
    freq_extrapolate(irr_params, general_params, rho_oveff_ratio_L1)
% freq_extrapolate
%
% Syntax:
%   [extrapolated_irr_params, rhof_veff_ratio_vector] = ...
%       freq_extrapolate(irr_params, sim_params, rho_oveff_ratio_L1)
%
% Description:
%   This function extrapolates the irregularity parameters (U, mu0, p1, p2) 
%   from a reference frequency (L1) to two other frequencies (L2 and L5). 
%   The scaling follows the conventions in [1]:
%       - Equation (17) for scaling the break wavenumber mu0.
%       - Equation (18) for scaling the ratio (rho_F / v_eff).
%       - Equation (19) for the piecewise definition of the universal 
%         turbulence strength parameter U.
%
% Inputs:
%   irr_params - Struct with fields:
%       .U   - Turbulence strength at L1
%       .mu0 - Normalized break wavenumber at L1
%       .p1  - Low-frequency spectral index
%       .p2  - High-frequency spectral index
%
%   sim_params  - Struct returned by get_sim_params(), containing 
%                        the field .gps_bands, a 1x3 vector with the frequencies 
%                        [L1, L2, L5].
%
%   rho_oveff_ratio_L1  - Ratio (rho_F / v_eff) known at L1 (scalar).
%
% Outputs:
%   extrapolated_irr_params - Struct with fields .L1, .L2 and .L5,
%       each containing .U, .mu0, .p1, and .p2. These represent the 
%       extrapolated parameters at the corresponding frequency.
%
%   rhof_veff_ratio_vector - 1x3 vector [rho_oveff_L1, rho_oveff_L2, rho_oveff_L5].
%
% Notes:
%   - This code is an adaptation from the function FreqExtrapolate.m available at:
%       https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/master/Libraries/GenScintFieldRealization/FreqExtrapolate.m
%
% Example:
%   % Suppose we have a single set of irregularity parameters:
%   irr_params.U   = 1.5;
%   irr_params.mu0 = 0.8;
%   irr_params.p1  = 2.0;
%   irr_params.p2  = 3.5;
%
%   % Get the general parameters (which includes gps_bands for L1, L2, L5)
%   sim_params = get_sim_params();
%
%   % Known ratio (rho_F / v_eff) at L1
%   rho_oveff_l1 = 1.0;
%
%   % Extrapolate
%   [extrap_params, rho_vec] = freq_extrapolate(irr_params, sim_params, rho_oveff_l1);
%
% Reference:
%   [1] Jiao, Yu, Rino, Charles, Morton, Yu (Jade), Carrano, Charles, 
%       "Scintillation Simulation on Equatorial GPS Signals for Dynamic 
%       Platforms," Proceedings of the 30th International Technical Meeting
%       of the Satellite Division of The Institute of Navigation (ION GNSS+
%       2017), Portland, Oregon, September 2017, pp. 1644-1657. 
%       https://doi.org/10.33012/2017.15258
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Extract frequencies (L1, L2, L5) from sim_params.
    freq_L1 = general_params.gps_bands(1);
    freq_L2 = general_params.gps_bands(2);
    freq_L5 = general_params.gps_bands(3);

    % Scale the reference ratio (rho_F / v_eff) for L2 and L5 [Eq. (13)].
    rho_oveff_ratio_L2 = rho_oveff_ratio_L1 * sqrt(freq_L1 / freq_L2);
    rho_oveff_ratio_L5 = rho_oveff_ratio_L1 * sqrt(freq_L1 / freq_L5);
    rhof_veff_ratio_vector = [rho_oveff_ratio_L1, rho_oveff_ratio_L2, rho_oveff_ratio_L5];

    % Extract parameters at L1.
    U_L1   = irr_params.U;
    mu0_L1 = irr_params.mu0;
    p1     = irr_params.p1;
    p2     = irr_params.p2;

    % Extrapolate from L1 to L2 and from L1 to L5.
    [U_L2, mu0_L2] = local_extrapolate(U_L1, mu0_L1, p1, p2, freq_L1, freq_L2);
    [U_L5, mu0_L5] = local_extrapolate(U_L1, mu0_L1, p1, p2, freq_L1, freq_L5);

    % Store results in the output struct.
    extrapolated_irr_params = struct('L1', [], 'L2', [], 'L5', []);

    extrapolated_irr_params.L1 = struct('U', U_L1,  ...
                                         'mu0', mu0_L1, ...
                                         'p1', p1, ...
                                         'p2', p2);
    extrapolated_irr_params.L2 = struct('U', U_L2,  ...
                                         'mu0', mu0_L2, ...
                                         'p1', p1, ...
                                         'p2', p2);
    extrapolated_irr_params.L5 = struct('U', U_L5,  ...
                                         'mu0', mu0_L5, ...
                                         'p1', p1, ...
                                         'p2', p2);
end