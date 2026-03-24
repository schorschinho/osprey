function indexMMLipBasisFunctionInBasis(obj)
%%  indexMMLipBasisFunctionInBasis(obj, input)
%   This method identifies MM/Lip entries in the basis set
%
%   USAGE:
%       obj.includeBasisFunctionInFit(obj)
%
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

    % Check which MM/LIP are available in the basis set
    % and match the input
    [ToInclude] = contains(obj.BasisSets.names,{'MM','Lip'});    % Get vector of logical indices
    ToInclude = obj.BasisSets.names(ToInclude);
    for rr = 1:length(ToInclude)
        idxToInclude = find(strcmp(ToInclude{rr}, obj.BasisSets.names));    % Get index of basis function to include
        obj.BasisSets.indexMMLipBasisFunction(obj.step+1,idxToInclude) = 1;                % Set index to 1
    end        
                                                                      % End loop over list of basis function names
end