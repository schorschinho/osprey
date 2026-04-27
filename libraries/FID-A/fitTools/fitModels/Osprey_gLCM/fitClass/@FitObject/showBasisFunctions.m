function tableOut = showBasisFunctions(obj,step)
%%  showBasisFunctions(obj)
%   This method returns a table showing the basis functions included in the
%   basis set, and whether they will be used in the fit.
%
%   USAGE:
%       obj.showBasisFunctions(step)
%
%   INPUTS:
%       step            = step to show    
%
%   OUTPUTS:
%       figure
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
%%  Diverge to default options if required 

    if nargin < 2
        step = obj.step;                                                % Set to last step
    end

%% Generate table

    tableOut = table(obj.BasisSets.names', obj.BasisSets.includeInFit(step,:)', ...
        'VariableNames', {'Basis function', 'Include in fit?'});
end