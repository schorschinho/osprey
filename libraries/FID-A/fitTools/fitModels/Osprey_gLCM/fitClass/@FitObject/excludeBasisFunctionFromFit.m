function excludeBasisFunctionFromFit(obj, input)
%%  excludeBasisFunctionFromFit(obj, input)
%   This method excludes basis functions according to the input
%
%   USAGE:
%       obj.excludeBasisFunctionFromFit(obj, input)
%
%   INPUTS:
%       input   = list of metabolites. %OPTIONS: 'all' or a list of names     
%       
%   OUTPUTS:
%       obj     = OspreyFitObj.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Generate index vector

    if ischar(input)                                                                % If is char create a 1x1 cell array
        input = {input};
    end
    
    if strcmpi(input, 'all')                                                        % If 'all' include all basis functions of the basis set
        input = obj.BasisSets.names;
    end
    
    for ii = length(input)                                                          % Loop over list of basis function names
        % Check which metabolites are available in the basis set
        % and match the input
        [metsToInclude, ~, ~] = intersect(obj.BasisSets.names, input, 'stable');    % Get vector metabolites in the basis set
        for rr = 1:length(metsToInclude)
            idxToExclude = find(strcmp(metsToInclude{rr}, obj.BasisSets.names));    % Get index of basis function to exclude
            obj.BasisSets.includeInFit(obj.step+1,idxToExclude) = 0;                % Set index to 0 
        end        
    end                                                                             % End loop over list of basis function names
end
        