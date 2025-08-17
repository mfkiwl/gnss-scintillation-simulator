function extrapolated_spectral_params = extrapolate_spectral_params(sim_params, freq)
%EXTRAPOLATE_SPECTRAL_PARAMS Extrapolates spectral parameters for scintillation simulation to a new frequency.
%
%   This function calculates the spectral parameters (p1, p2, U, mu0) for a
%   given target frequency based on a set of reference parameters. The
%   extrapolation follows the methodology described for simulating GPS
%   scintillation signals on dynamic platforms, particularly the frequency
%   scaling laws for the universal scattering strength (U) and the
%   normalized break wavenumber (mu0).
%
%   The function first extrapolates the normalized break wavenumber (mu0)
%   and then uses it, along with the reference parameters, to determine the
%   correct case for extrapolating the universal scattering strength (U).
%
%   Syntax:
%   extrapolated_spectral_params = extrapolate_spectral_params(sim_params, freq)
%
%   Inputs:
%   sim_params - A structure containing the simulation parameters. It must
%                include reference spectral parameters for a given severity
%                level (e.g., 'strong', 'moderate'). The required fields are:
%                .const.spectral.freq_ref.value - Reference frequency (Hz).
%                .severity - A string indicating the scintillation severity.
%                .const.spectral.(severity).p1 - Spectral index 1.
%                .const.spectral.(severity).p2 - Spectral index 2.
%                .const.spectral.(severity).U_ref - Reference universal
%                                                   scattering strength.
%                .const.spectral.(severity).mu0_ref - Reference normalized
%                                                     break wavenumber.
%
%   freq       - The target frequency (Hz) to which the spectral
%                parameters will be extrapolated.
%
%   Outputs:
%   extrapolated_spectral_params - A structure containing the extrapolated
%                                  spectral parameters for the target frequency.
%   .p1   - Spectral index 1 [1].
%   .p2   - Spectral index 2 [1].
%   .U    - Extrapolated universal scattering strength.
%   .mu0  - Extrapolated normalized break wavenumber.
%
%   Refs:
%   [1] Jiao, Y., Rino, C., Morton, Y. J., & Carrano, C. (2017).
%   Scintillation Simulation on Equatorial GPS Signals for Dynamic Platforms.
%   [2] Xu, Dongyang, Y.T. Jade Morton, Charles L. Rino, Charles S. Carrano, and Yu Jiao.
%   “A Two-Parameter Multifrequency GPS Signal Simulator for Strong Equatorial 
%   Ionospheric Scintillation: Modeling and Parameter Characterization.” NAVIGATION 67, 
%   no. 1 (2020): 181–95. https://doi.org/10.1002/navi.350.
%
%% Initialization
freq_ref = sim_params.const.spectral.freq_ref.value;
severity = sim_params.severity;

p1       = sim_params.const.spectral.(severity).p1;
p2       = sim_params.const.spectral.(severity).p2;
U_ref    = sim_params.const.spectral.(severity).U_ref ;
mu0_ref  = sim_params.const.spectral.(severity).mu0_ref;

if freq == freq_ref
    extrapolated_spectral_params.p1  = p1;
    extrapolated_spectral_params.p2  = p2;
    extrapolated_spectral_params.U   = U_ref;
    extrapolated_spectral_params.mu0 = mu0_ref;
    % early return
    return
end

%% Extrapolate the spectral parameter μ₀: the normalized break wavenumber

% Extrapolate mu0 from freq_ref to freq based on Eq. (15)
mu0 = mu0_ref * sqrt(freq_ref / freq);

%% Extrapolate the spectral parameter U

% Exponents for (freq_ref/freq) based on p1.
exponent_p1 = 0.5 * p1 + 1.5;

% Determine the cases for U(freq) computation.
cond_ref = (mu0_ref >= 1);
cond_freq = (mu0 >= 1);

% SEE: [1, Eq. (17)]
% BUG: For some scaterring regime and GNSS frequency combination, the simulation fails. I am not sure that [1] is the best reference. I didn't read it thoroughly, but [2, Eq. (14)] seems to extrapolate U differently. We should understand the theoretical difference between both approches and use [2] instead of [1] if it is suitable.
if cond_ref && cond_freq
    % Case 1: mu0(freq_ref) >= 1, mu0(freq) >= 1
    U = U_ref * (freq_ref/freq)^exponent_p1;
elseif ~cond_ref && cond_freq
    % Case 2: mu0(freq_ref) < 1, mu0(freq) >= 1
    U = U_ref / mu0_ref^(p2 - p1) * (freq_ref/freq)^exponent_p1;
else
    % Case 3: mu0(freq_ref) < 1, mu0(freq) < 1
    assert((mu0_ref < 1) && (mu0 < 1), ['The scattering regime (%s) and signal frequency (%d MHz) combination has led to\n' ...
        'improper conditions for scattering strength (U) extrapolation. That is probably a theoterical\n' ...
        'limitation rather than a code error. Please refer to "Scintillation Simulation on Equatorial\n' ...
        'GPS Signals for Dynamic Platforms", from Yu Jiao et. al., and check the source code to check\n' ...
        'whether there is a solution for that. If not, you should avoid using this frequency and\n' ...
        'scattering regime combination.'], severity, freq);
    U = U_ref * (mu0 / mu0_ref)^(p2 - p1) * (freq_ref/freq)^exponent_p1;
end

%% spectral parameters output
extrapolated_spectral_params.p1 = p1;
extrapolated_spectral_params.p2 = p2;
extrapolated_spectral_params.U = U;
extrapolated_spectral_params.mu0 = mu0;

end
