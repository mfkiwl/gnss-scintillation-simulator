function rad = dms_str_2_rad(dmsStr)
% dms_str_2_rad
%
% Syntax:
%   rad = dms_str_2_rad(dmsStr)
%
% Description:
%   This function converts a coordinate given as a string in DMS 
%   (degree-minute-second) format into radians. It parses the input string,
%   converts the degrees, minutes, and seconds to a decimal degree value,
%   and then converts that value to radians.
%
% Inputs:
%   dmsStr - A string representing the coordinate in DMS format. The format
%            must be 'deg°min''sec"Direction', where Direction is one of
%            'N', 'S', 'E', or 'W'.
%
% Outputs:
%   rad - The coordinate in radians.
%
% Notes:
%   - The function expects the DMS string to match the specified format
%     exactly, including the degree (°), minute ('), and second (") symbols.
%   - An error is thrown if the input string cannot be parsed.
%
% Example:
%   % Convert '45°30''15"N' to radians:
%   radians = dms_str_2_rad('45°30''15"N');
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Parse the DMS string using regular expressions.
    tokens = regexp(dmsStr, '(\d+)°(\d+)''([\d\.]+)"?([NSWE])', 'tokens');
    if isempty(tokens)
        error('Could not parse the DMS string.');
    end
    parts = tokens{1};
    
    % Convert degrees, minutes, and seconds to numbers.
    degVal = str2double(parts{1});
    minVal = str2double(parts{2});
    secVal = str2double(parts{3});
    
    % Determine sign: negative for South or West.
    if any(strcmpi(parts{4}, {'S','W'}))
        signVal = -1;
    else
        signVal = 1;
    end
    
    % Convert to decimal degrees.
    decimalDegrees = signVal * (degVal + minVal/60 + secVal/3600);
    
    % Convert decimal degrees to radians.
    rad = deg2rad(decimalDegrees);
end
