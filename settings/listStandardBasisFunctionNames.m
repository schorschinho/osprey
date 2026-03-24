% listStandardBasisFunctionNames.m
% Georg Oeltzschner, Johns Hopkins University 2022
%
% USAGE:
% allMets = listStandardBasisFunctionNames
%
% DESCRIPTION:
% This function returns a cell array of standard metabolite and
% macromolecular basis function names. It will look up the list of VALID
% names via 'listValidBasisFunctionNames.m'. For metabolites with multiple
% valid names, it will only return the first (standard) name.
%
% OUTPUTS:
% allMets = Cell array with a list of standard metabolite/MM names according
%           to Osprey naming convention.
%
% INPUTS
% type    = String (can be 'mets' or 'mm')

function listOfStandardNames = listStandardBasisFunctionNames(type)

% Get the list of valid names
listOfStandardNames = listValidBasisFunctionNames(type);

% Loop over elements and remove all non-default names
for rr = 1:length(listOfStandardNames)
    if iscell(listOfStandardNames{rr})
        % Pick first element (which is the default name)
        listOfStandardNames{rr} = listOfStandardNames{rr}{1};
    else
        % No need to do anything
    end
    
end

end
