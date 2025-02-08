function [propagated_scint_field,norm_phase_sdf,detrended_phase_realization,mu,doppler_frequency] = get_scintillation_time_series(gen_params, irr_params, rhof_veff_ratio, seed)
% get_scintillation_time_series
%
% Syntax:
%   propagated_scint_field = get_scintillation_time_series(gen_params, ...
%                                                          irr_params, ...
%                                                          rhof_veff_ratio, ...
%                                                          seed)
%
% Description:
%   Generates a time-domain scintillation field by performing the following steps:
%   1) Determines a suitable FFT length (nfft) based on the simulation time and 
%      sampling interval (dt) via the helper function nicefftnum.m.
%   2) Computes a Doppler frequency axis over the range [-Nyquist, +Nyquist).
%   3) Constructs the spatial frequency axis (mu) by scaling the Doppler 
%      axis with 2*pi*rhof_veff_ratio.
%   4) Obtains the normalized phase spectral density function (norm_phase_sdf) 
%      from get_norm_phase_sdf, using the irregularity parameters (irr_params) 
%      and mu.
%   5) Calls get_phase_realization to generate a random phase realization 
%      based on norm_phase_sdf.
%   6) Propagates the resulting phase-perturbed wavefield by calling 
%      get_propagated_field, yielding the final time-domain scintillation 
%      field.
%
% Inputs:
%   gen_params - Struct containing simulation parameters, including:
%       .simulation_time : Total simulation time (seconds)
%       .dt              : Time step (seconds)
%
%   irr_params - Struct with irregularity parameters required by 
%                get_norm_phase_sdf, typically containing:
%       .U               : Turbulence strength
%       .mu0             : Break wavenumber
%       .p1, .p2         : Spectral indices
%
%   rhof_veff_ratio    - Scalar representing (rho_F / v_eff), used to scale 
%                        the Doppler frequency axis into a spatial/temporal 
%                        frequency axis mu.
%
%   seed               - Integer seed for the random number generator used 
%                        by get_phase_realization.
%
% Outputs:
%   propagated_scint_field - Complex time-domain scintillation field after 
%                            applying phase perturbations and parabolic wave 
%                            propagation.
%
% Dependencies:
%   - nicefftnum(sim_time_ratio)
%       Calculates an FFT size based on the ratio of simulation_time/dt.
%   - get_norm_phase_sdf(mu, irr_params)
%       Computes a normalized phase spectral density function given mu and 
%       irregularity parameters.
%   - get_phase_realization(norm_phase_sdf, D_mu, seed)
%       Generates a random (possibly complex) phase realization using the 
%       normalized phase spectral density. Requires a parameter D_mu not shown 
%       in this snippet.
%   - get_propagated_field(mu, detrended_phase_realization)
%       Applies a parabolic propagation factor in the frequency domain and 
%       returns the propagated scintillation field.
%
% Notes:
%   - The development of this function was inspired by the code available at:
%     https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param/blob/master/Libraries/GenScintFieldRealization/GenScintFieldRealization.m
%
% Example:
%   % Assuming gen_params, irr_params, and D_mu are already defined or loaded:
%   gen_params.simulation_time = 60;  % seconds
%   gen_params.dt = 0.01;            % time step
%   irr_params.U   = 1.5;             % turbulence strength
%   irr_params.mu0 = 0.8;             % break wavenumber
%   irr_params.p1  = 2.0;             % spectral index (low freq)
%   irr_params.p2  = 3.5;             % spectral index (high freq)
%   ratio         = 0.5;             % example (rho_F / v_eff)
%   seed_val      = 12345;           % random seed
%
%   scint_field = get_scintillation_time_series(gen_params, ...
%                                               irr_params, ...
%                                               ratio, ...
%                                               seed_val);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    nfft = nicefftnum(gen_params.simulation_time / gen_params.dt);
    doppler_frequency = (-nfft/2 : nfft/2-1) / (nfft * gen_params.dt); 
    mu = 2 * pi * doppler_frequency * rhof_veff_ratio;
    D_mu = mu(2) - mu(1);

    % Obtain the normalized phase spectral density function.
    norm_phase_sdf = get_norm_phase_sdf(mu, irr_params);

    % Generate the random phase realization (note: 'D_mu' must be defined externally).
    detrended_phase_realization = get_phase_realization(norm_phase_sdf, D_mu, nfft, seed);

    % Propagate the scintillation field.
    propagated_scint_field = get_propagated_field(mu, detrended_phase_realization);

end
