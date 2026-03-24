% fit_sortBasisSet.m
% Helge Zoellner, Johns Hopkins University 2020.
%
% USAGE:
% basisSetOut = fit_sortBasisSet(basisSetIn)
% 
% DESCRIPTION:
% This function sorts an external basis set according to the conventions
% used in Osprey
% 
% OUTPUTS:
% basisSetOut = FID-A basis set container that only contains the basis
%               functions specified in metabList
%
% INPUTS:
% basisSetIn  = FID-A basis set container (loaded with io_LCMBasis).

function basisSetOut = fit_sortBasisSet(basisSetIn)

% Duplicate the input basis set
basisSetOut = basisSetIn;

% Retrieve the list of valid metabolite names
allMets = listValidBasisFunctionNames('mets');

% Loop over all entries in the list. If there is a match with the list of
% names in the input basis set, add the respective basis function to the
% output basis set.
actBasisFnct = 1;
for kk = 1 : length(allMets)
    name = allMets{kk};
    % If there is more than one synonym, search for all of them
    idx          = find(ismember(basisSetIn.name,name));
    if ~isempty(idx)
        % If there is a match, name accordingly to the list of valid
        % metabolite names
        if iscell(name)
            % If the metabolite has several synonyms, name it according
            % to the first entry (which is the Osprey default).
            basisSetOut.name{actBasisFnct}   = name{1};
        else
            basisSetOut.name{actBasisFnct}   = name;
        end
        basisSetOut.fids(:,actBasisFnct,:)   = basisSetIn.fids(:,idx,:);
        basisSetOut.specs(:,actBasisFnct,:)  = basisSetIn.specs(:,idx,:);
        actBasisFnct = actBasisFnct + 1;
    end
end
basisSetOut.nMets = actBasisFnct-1;

% Do the same for MMs
allMMs = listValidBasisFunctionNames('mm');
for kk = 1 : length(allMMs)
    name = allMMs{kk};
    % If there is more than one synonym, search for all of them
    idx          = find(ismember(basisSetIn.name,name));
    if ~isempty(idx)
        % If there is a match, name accordingly to the list of valid
        % metabolite names
        if iscell(name)
            % If the metabolite has several synonyms, name it according
            % to the first entry (which is the Osprey default).
            basisSetOut.name{actBasisFnct}   = name{1};
        else
            basisSetOut.name{actBasisFnct}   = name;
        end
        basisSetOut.fids(:,actBasisFnct,:)   = basisSetIn.fids(:,idx,:);
        basisSetOut.specs(:,actBasisFnct,:)  = basisSetIn.specs(:,idx,:);
        actBasisFnct = actBasisFnct + 1;
    end
end
basisSetOut.nMM = actBasisFnct - basisSetOut.nMets-1;

% Since we duplicated the input basis set in the beginning and just kept
% the valid metabolites that we could find in the basis set, we need to
% remove any 'surplus' ones at the end.
try
    basisSetOut.name(actBasisFnct:end) = [];
    basisSetOut.fids(:,actBasisFnct:end,:)   = [];
    basisSetOut.specs(:,actBasisFnct:end,:)  = [];
catch
end

% Write new size of the basis set.
basisSetOut.sz = size(basisSetOut.fids);

end