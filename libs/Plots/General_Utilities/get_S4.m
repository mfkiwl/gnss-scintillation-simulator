% get_S4
% Computes the S4 scintillation index from the signal intensity (SI).
%
% Syntax:
%   S4 = get_S4(SI)
%
% Description:
%   This function calculates the S4 scintillation index, a widely used metric 
%   to quantify ionospheric scintillation effects on GNSS signals. The S4 index 
%   is defined as the normalized standard deviation of the received signal 
%   intensity (SI), capturing the magnitude of intensity fluctuations.
%
% Inputs:
%   SI - A vector of signal intensity values representing the received signal 
%        power variations over time.
%
% Outputs:
%   S4 - The computed S4 scintillation index, given by:
%        S4 = sqrt((mean(SI.^2) - mean(SI)^2) / (mean(SI)^2))
%
% Notes:
%   - A higher S4 value indicates stronger intensity scintillation effects.
%   - The function assumes that SI is a real-valued vector representing 
%     the signal intensity fluctuations in a GNSS receiver.
%
% References:
%   [1] Vasylyev D, Béniguel Y, Volker W, Kriegel M & Berdermann J, et al. 
%     2022. Modeling of ionospheric scintillation. J. Space Weather 
%     Space Clim. 12, 22. https://doi.org/10.1051/swsc/2022016.
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 ORCID: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com


function S4 = get_S4(SI)
    % Equation (1) of [1].
    S4 = sqrt((mean(SI.^2) - mean(SI)^2)/(mean(SI)^2));
end