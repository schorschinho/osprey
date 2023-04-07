% fit_readLCMBasisSetMetabs.m
% Georg Oeltzschner, Johns Hopkins University 2023.
%
% USAGE:
% metabList = fit_readLCMBasisSetMetabs(filename);
%
% DESCRIPTION:
% Reads an LCModel basis set and returns a cell array with all metabolite
% names in the basis set.
%
% INPUTS:
% filename   = filename of LCModel .basis file.
%
% OUTPUTS:
% metabList = cell array with all metabolite names

function metabList = fit_readLCMBasisSetMetabs(filename)

% Begin to read the basis set file.
fid     = fopen(filename);
line    = fgets(fid);

%look for lines containing the variable METABO
metabList = {};
while ~feof(fid)
    
    % Read the key-value pair
    [key, value] = readKeyValuePair(line);
    if ~isempty(key)
        if strcmpi(key, 'METABO')
            % If a metabolite is found, run it against the metabolite list
            value = strrep(value, ',', '');
            value = strrep(value, '''', '');
            metabList{end+1} = value;
        end
    end
    
    line = fgets(fid);

end

fclose(fid);

end

% Helper function
function [key, value] = readKeyValuePair(line)
    % Use a regular expression to extract the key and value
    expression = '([^\s=]+)\s*=\s*([^\s=]+)';
    [tokens, match] = regexp(line, expression, 'tokens', 'match');
    
    % If a match was found, assign the key and value to the output variables
    if ~isempty(match)
        key = tokens{1}{1};
        value = tokens{1}{2};
    else
        % If no match was found, set the key and value to empty strings
        key = '';
        value = '';
    end
end
