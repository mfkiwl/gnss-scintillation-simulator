function [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2)
% local_extrapolate
%
% Syntax:
%   [U_f2, mu0_f2] = local_extrapolate(U_f1, mu0_f1, p1, p2, f1, f2)
%
% Description:
%   Applies [1, Eq. 17] to scale mu0 between frequencies, and applies 
%   the piecewise expression of [1, Eq. 19] to determine the new universal 
%   turbulence strength U. The specific case depends on whether mu0 >= 1 
%   at the reference frequency (f1) or at the target frequency (f2).
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

    % Exponents for (f1/f2) based on p1.
    exponent_p1 = 0.5 * p1 + 1.5;

    % Scale mu0 from f1 to f2.
    mu0_f2 = mu0_f1 * sqrt(f1 / f2);

    % Determine the cases for U_f2 computation.
    cond_f1 = (mu0_f1 >= 1);
    cond_f2 = (mu0_f2 >= 1);

    if cond_f1 && cond_f2
        % Case 1: mu0(f1) >= 1, mu0(f2) >= 1
        U_f2 = U_f1 * (f1/f2)^exponent_p1;
    elseif ~cond_f1 && cond_f2
        % Case 2: mu0(f1) < 1, mu0(f2) >= 1
        U_f2 = U_f1 / mu0_f1^(p2 - p1) * (f1/f2)^exponent_p1;
    else
        % Case 3: mu0(f1) < 1, mu0(f2) < 1
        U_f2 = U_f1 * (mu0_f2 / mu0_f1)^(p2 - p1) * (f1/f2)^exponent_p1;
    end
end
