function rhof_veff_ratio = extrapolate_scaling_param(sim_params, freq, rhof_veff_ratio_ref)
%freq_rhof_veff Extrapolate scaling parameter to a new frequency.
%
%   rhof_veff_ratio = freq_rhof_veff(sim_params, freq, rhof_veff_ratio_ref)
%   scales the reference ratio (rho_F/v_eff) at the model's reference
%   frequency to a target frequency using [1, Eq. (13)].
%
% Inputs:
%   sim_params              - (required, struct) Simulation parameters with fields:
%                              cte.spectral.freq_ref.value (reference freq in Hz).
%   freq                    - (required, Hz, scalar) Target frequency for scaling.
%   rhof_veff_ratio_ref     - (required, scalar) Reference rho_F/v_eff ratio.
%
% Output:
%   rhof_veff_ratio         - (scalar) Scaled rho_F/v_eff ratio at target frequency.
%
% Reference:
%  [1] Jiao, Yu, Rino, Charles, Morton, Yu (Jade), Carrano, Charles,
%  "Scintillation Simulation on Equatorial GPS Signals for Dynamic
%  Platforms," Proceedings of the 30th International Technical Meeting
%  of the Satellite Division of The Institute of Navigation (ION GNSS+
%  2017), Portland, Oregon, September 2017, pp. 1644-1657.
%  https://doi.org/10.33012/2017.15258
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com
%
%   Rubem Vasconcelos Pacelli
%   ORCID: https://orcid.org/0000-0001-5933-8565
%   Email: rubem.engenharia@gmail.com

%% Initialization
freq_ref = sim_params.cte.spectral.freq_ref.value;

%% Extrapolate the scaling parameter
% Scale the reference ratio (rho_F / v_eff) for L2 and L5 [1, Eq. (13)].
rhof_veff_ratio = rhof_veff_ratio_ref * sqrt(freq_ref / freq);
end