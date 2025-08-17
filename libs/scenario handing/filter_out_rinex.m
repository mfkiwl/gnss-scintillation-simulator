function rinex = filter_out_rinex(log, rinex, prn, constellation, frequency)
%FILTER_OUT_RINEX Summary of this function goes here
%   Detailed explanation goes here

%% Filter out undesired constellations
% if the user passed no constellation, they must decide on that at runtime
% interactive constellation selection
if isempty(char(constellation))
    fprintf(['You have not passed any constellation. ' ...
        'For this RINEX file, the following are availables:\n\n']);

    rinex_fieldnames = fieldnames(rinex);
    for i = 1:numel(rinex_fieldnames)
        fprintf('%d) %s\n', i, rinex_fieldnames{i});
    end
    select = input(['\nPass a comma-separated number sequence ' ...
        '(no spaces in between) to select which constellation(s) ' ...
        'you want\n'], "s");
    
    % get the selected indices
    select = cellfun(@(x) double(string(x)), strsplit(select, ','));
    % create a string array with the selected constellation(s)
    constellation = cellfun(@(x) lower(string(x)), ...
        rinex_fieldnames(select));
end

% ensure that constellation is a string and not a char array
constellation = string(constellation);

rinex_fieldnames = cellfun(@string, fieldnames(rinex));
fields_to_remove = ~ismember(lower(rinex_fieldnames), constellation);
if all(fields_to_remove)
    log.error('There is no constellation(s) %s in this RINEX file.', ...
        strjoin(constellation, ", "));
end
% remove undesired constellation fields
rinex = rmfield(rinex, rinex_fieldnames(fields_to_remove));

%% Filter out undesired PRNs
if ~isempty(char(prn))
    % TODO: filter out PRNs
end

%% Filter out undesired frequencies


if ~isempty(char(frequency))
    % TODO: filter out frequencies
end

end