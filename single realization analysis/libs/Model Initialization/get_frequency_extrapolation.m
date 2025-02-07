function [extrapolated_irregularity_parameters, extrapolated_rhof_veff_ratio] = ...
    freq_extrapolate_all(irregularity_params, rho_oveff_f1, f1, f2)
% FREQ_EXTRAPOLATE_ALL  Applies frequency-extrapolation (Eqs. (12), (13), (14))
%                       to all conditions (Severe, Moderate, Mild).
%
% Syntax:
%   [extrapolated_irregularity_parameters, extrapolated_rhof_veff_ratio] ...
%       = freq_extrapolate_all(irregularity_params, rho_oveff_f1, f1, f2)
%
% Description:
%   This function automatically iterates over the 'Severe', 'Moderate', and 'Mild'
%   fields of the input `irregularity_params` struct, applying the scaling
%   relationships (Eqs. (12), (13), (14) in the paper) for U and mu0 from 
%   frequency f1 to frequency f2. The ratio rho_F / v_eff is assumed to be
%   the same across conditions, so only a single new value for rho_F / v_eff
%   is computed.
%
% Inputs:
%   irregularity_params  - Struct with subfields .Severe, .Moderate, .Mild, each
%                          containing [U, mu0, p1, p2]. Obtained from 
%                          get_irregularity_parameters().
%
%   rho_oveff_f1         - The ratio (rho_F / v_eff) at frequency f1 (scalar).
%
%   f1, f2               - Original and target frequencies (Hz).
%
% Outputs:
%   extrapolated_irregularity_parameters - Struct with the same fields (.Severe,
%                                          .Moderate, .Mild), but with U and mu0
%                                          scaled to frequency f2. (p1, p2 remain 
%                                          unchanged.)
%
%   extrapolated_rhof_veff_ratio         - Single scalar with the new rho_F / v_eff 
%                                          at f2, computed via Eq. (13).
%
% Equations (refer to the paper):
%   (12) mu0(f_s^i, f_r) = mu0(f_r) * sqrt(f_r / f_s)
%   (13) rho_F / v_eff(f_s^i, f_r) = rho_F / v_eff(f_r) * sqrt(f_r / f_s)
%   (14) Piecewise definitions for U(f_s^i, f_r) depending on mu0 >= 1 or < 1.
%
% Example:
%   irr_params = get_irregularity_parameters(); % default param set
%   f1 = 1.57542e9;    % e.g., GPS L1 freq
%   f2 = 1.22760e9;    % e.g., GPS L2 freq
%   rho_oveff_L1 = 1.0;
%
%   [extrap_params, new_rho] = freq_extrapolate_all(irr_params, rho_oveff_L1, f1, f2);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Pre-allocate the output struct with the same fields as the input.
    extrapolated_irregularity_parameters = struct('Severe', [], 'Moderate', [], 'Mild', []);

    % Since rho_F / v_eff does not depend on the particular irregularity parameters,
    % we compute it once using Eq. (13).
    extrapolated_rhof_veff_ratio = rho_oveff_f1 * sqrt(f1 / f2);

    % List of condition fields to iterate.
    condition_list = {'Severe', 'Moderate', 'Mild'};

    % Loop through each condition in the irregularity_params struct.
    for condition = condition_list
        cond_str = condition{1};  % Extract the actual string ('Severe', 'Moderate', or 'Mild')

        % Extract U, mu0, p1, p2 from the chosen sub-struct.
        U   = irregularity_params.(cond_str).U;
        mu0 = irregularity_params.(cond_str).mu0;
        p1  = irregularity_params.(cond_str).p1;
        p2  = irregularity_params.(cond_str).p2;

        %% (1)  Apply Eq. (14) logic at f1 => intermediate Cpp
        if mu0 >= 1
            Cpp = U;
        else
            Cpp = U / mu0^(p2 - p1);
        end

        %% (2)  Scale Cpp from f1 to f2 using Joy's exponent: (f1/f2)^(0.5*p1 + 1.5)
        Cpp_f2 = Cpp * (f1 / f2)^(0.5 * p1 + 1.5);

        %% (3)  mu0(f2) via Eq. (12).
        mu0_f2 = mu0 * sqrt(f1 / f2);

        %% (4)  U_f2 via Eq. (14) logic at f2.
        if mu0_f2 >= 1
            U_f2 = Cpp_f2;
        else
            U_f2 = Cpp_f2 * (mu0_f2)^(p2 - p1);
        end

        % Create the updated sub-struct for this condition.
        extrapolated_irregularity_parameters.(cond_str) = struct( ...
            'U',   U_f2,    ...
            'mu0', mu0_f2,  ...
            'p1',  p1,      ...
            'p2',  p2       ...
        );
    end
end
