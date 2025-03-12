function rad = dms_str_2_rad(dmsStr)
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
