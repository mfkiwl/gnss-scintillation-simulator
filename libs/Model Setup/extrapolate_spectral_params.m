function extrapolated_spectral_params = extrapolate_spectral_params(sim_params, freq)
%EXTRAPOLATE_SPECTRAL_PARAMS Summary of this function goes here
%   Detailed explanation goes here
%% Initialization
freq_ref = sim_params.cte.spectral.freq_ref.value;
severity = sim_params.severity;

p1       = sim_params.cte.spectral.(severity).p1;
p2       = sim_params.cte.spectral.(severity).p2;
U_ref    = sim_params.cte.spectral.(severity).U_ref ;
mu0_ref  = sim_params.cte.spectral.(severity).mu0_ref;

if freq == freq_ref
    extrapolated_spectral_params.p1  = p1;
    extrapolated_spectral_params.p2  = p2;
    extrapolated_spectral_params.U   = U_ref;
    extrapolated_spectral_params.mu0 = mu0_ref;
    % early return
    return
end

%% Extrapolate the spectral parameter μ₀: the normalized TODO:...

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
        'going on and assuming that the extrapolation step is right.']);
    U = U_ref * (mu0 / mu0_ref)^(p2 - p1) * (freq_ref/freq)^exponent_p1;
end

%% spectral parameters output
extrapolated_spectral_params.p1 = p1;
extrapolated_spectral_params.p2 = p2;
extrapolated_spectral_params.U = U;
extrapolated_spectral_params.mu0 = mu0;

end

