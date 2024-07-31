function [basisSet] = scrubMMLipFromBasis(basisSet, verbose)

% Georg Oeltzschner, Johns Hopkins University 2023
% Removes MM and lipid basis functions from input basis set.

if nargin < 2
    verbose = 0;
end

% Go through list of valid MM names
listOfValidNames = listValidBasisFunctionNames('mm');
% Expand list (because it contains some 'multiple' entries that describe
% commonly used name for the same MM component
for ll = 1:length(listOfValidNames)
    if ~ischar(listOfValidNames{ll})
        for mm = 1:length(listOfValidNames{ll})
            listOfValidNames{end+1} = listOfValidNames{ll}{mm};
        end
        listOfValidNames(ll) = [];
    end
end

% Determine which entries of the list of valid MM names are matched in the
% basis set name array
whichIndicesAreTheMMs = ismember(basisSet.name, listOfValidNames);

% Report
if verbose
    fprintf('scrubMMLipFromBasis removes the following %i MM/Lip basis functions from the basis set: \n', sum(whichIndicesAreTheMMs));
    fprintf('%s\n', basisSet.name{whichIndicesAreTheMMs});
end

% Retain only basis functions that are MMs
basisSet.name   = basisSet.name(~whichIndicesAreTheMMs);
basisSet.fids   = basisSet.fids(:,~whichIndicesAreTheMMs,:);
basisSet.specs  = basisSet.specs(:,~whichIndicesAreTheMMs,:);

% Save metadata
basisSet.sz     = size(basisSet.fids);
basisSet.nMets  = sum(~whichIndicesAreTheMMs);
basisSet.nMM    = 0;

end