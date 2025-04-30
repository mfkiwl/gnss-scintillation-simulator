function rhof_veff_ratio = freq_rhof_veff(sim_params, freq, rhof_veff_ratio_ref)
% freq_extrapolate
%
% Syntax:
%   [extrapolated_irr_params, rhof_veff_ratio_vector] = ...
%       freq_extrapolate(irr_params, sim_params, rho_oveff_ratio_freq_ref)
%
% Description:
%   This function extrapolates the irregularity parameters (U, mu0, p1, p2)
%   from a reference frequency (freq_ref) to two other frequencies (L2 and L5).
%   The scaling follows the conventions in [1]:
%       - Equation (17) for scaling the break wavenumber mu0.
%       - Equation (18) for scaling the ratio (rho_F / v_eff).
%       - Equation (19) for the piecewise definition of the universal
%         turbulence strength parameter U.
%
% Inputs:
%   irr_params - Struct with fields:
%       .U   - Turbulence strength at freq_ref
%       .mu0 - Normalized break wavenumber at freq_ref
%       .p1  - Low-frequency spectral index
%       .p2  - High-frequency spectral index
%
%   sim_params  - Struct returned by get_sim_params(), containing
%                        the field .gps_bands, a 1x3 vector with the frequencies
%                        [freq_ref, L2, L5].
%
%   rho_oveff_ratio_freq_ref  - Ratio (rho_F / v_eff) known at freq_ref (scalar).
%
% Outputs:
%   extrapolated_irr_params - Struct with fields .freq_ref, .L2 and .L5,
%       each containing .U, .mu0, .p1, and .p2. These represent the
%       extrapolated parameters at the corresponding frequency.
%
%   rhof_veff_ratio_vector - 1x3 vector [rho_oveff_freq_ref, rho_oveff_L2, rho_oveff_L5].
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
%   % Get the general parameters (which includes gps_bands for freq_ref, L2, L5)
%   sim_params = get_sim_params();
%
%   % Known ratio (rho_F / v_eff) at freq_ref
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

%% Initialization
freq_ref = sim_params.cte.spectral.freq_ref.value;

%% Extrapolate the scaling parameter
% Scale the reference ratio (rho_F / v_eff) for L2 and L5 [Eq. (13)].
rhof_veff_ratio = rhof_veff_ratio_ref * sqrt(freq_ref / freq);
end