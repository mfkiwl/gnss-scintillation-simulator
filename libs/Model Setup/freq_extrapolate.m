function [rho_veff_ratio, extrapolated_spectral_params] = ...
    freq_extrapolate(spectral_params_ref, rho_veff_ratio_ref, freq_ref, freq)
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

%% freq_ref == freq -> No extrapolation
if freq == freq_ref
    rho_veff_ratio = rho_veff_ratio_ref;

    extrapolated_spectral_params.p1  = spectral_params_ref.p1;
    extrapolated_spectral_params.p2  = spectral_params_ref.p2;
    extrapolated_spectral_params.U   = spectral_params_ref.U_ref;
    extrapolated_spectral_params.mu0 = spectral_params_ref.mu0_ref;
    % early return
    return
end

%% Initialization
% Extract parameters at the reference frequency
U_ref   = spectral_params_ref.U_ref;
mu0_ref = spectral_params_ref.mu0_ref;
p1     = spectral_params_ref.p1;
p2     = spectral_params_ref.p2;

%% Extrapolate the spectral parameter μ₀

% Extrapolate mu0 from freq_ref to freq
mu0 = mu0_ref * sqrt(freq_ref / freq);

%% Extrapolate the spectral parameter U

% Exponents for (freq_ref/freq) based on p1.
exponent_p1 = 0.5 * p1 + 1.5;

% Determine the cases for U(freq) computation.
cond_ref = (mu0_ref >= 1);
cond_freq = (mu0 >= 1);

if cond_ref && cond_freq
    % Case 1: mu0(freq_ref) >= 1, mu0(freq) >= 1
    U = U_ref * (freq_ref/freq)^exponent_p1;
elseif ~cond_ref && cond_freq
    % Case 2: mu0(freq_ref) < 1, mu0(freq) >= 1
    U = U_ref / mu0_ref^(p2 - p1) * (freq_ref/freq)^exponent_p1;
else
    % Case 3: mu0(freq_ref) < 1, mu0(freq) < 1
    assert((mu0_ref < 1) && (mu0 < 1), ['Unexpected values of μ₀ for both reference ' ...
        'and target frequencies for extrapolation of U. You should check\n' ...
        'why the source code is wrongly reaching at this point before\n' ...
        'going on and assuming that the extrapolation spet is right.']);
    U = U_ref * (mu0 / mu0_ref)^(p2 - p1) * (freq_ref/freq)^exponent_p1;
end

%% spectral parameters output
extrapolated_spectral_params.p1 = spectral_params_ref.p1;
extrapolated_spectral_params.p2 = spectral_params_ref.p2;
extrapolated_spectral_params.U = U;
extrapolated_spectral_params.mu0 = mu0;

%% Extrapolate the scaling parameter
% Scale the reference ratio (rho_F / v_eff) for L2 and L5 [Eq. (13)].
rho_veff_ratio = rho_veff_ratio_ref * sqrt(freq_ref / freq);
end