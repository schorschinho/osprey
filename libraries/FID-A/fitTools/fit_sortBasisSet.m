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

% Save all available metabolite names in a cell
all_mets = {'Ala','Asc','Asp','bHB','bHG','Cit','Cr','CrCH2','EA','EtOH','GABA','GPC','GSH','Glc','Gln' ...
    ,'Glu','Gly','H2O','Ins','Lac','NAA','NAAG','PCh','PCr','PE','Phenyl' ...
    ,'Scyllo','Ser','Tau','Tyros'};

% Duplicate the input basis set
basisSetOut = basisSetIn;

actBasisFnct = 1;
for kk = 1 : length(all_mets)
    name = all_mets{kk};
    idx          = find(strcmp(basisSetIn.name,name));
    if ~isempty(idx)
        basisSetOut.name{actBasisFnct} = basisSetIn.name{idx};
        basisSetOut.fids(:,actBasisFnct,:)   = basisSetIn.fids(:,idx,:);
        basisSetOut.specs(:,actBasisFnct,:)  = basisSetIn.specs(:,idx,:);
        actBasisFnct = actBasisFnct + 1;
    end
end
basisSetOut.nMets = actBasisFnct-1;

all_MMs = {'MM09','MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20'};
for kk = 1 : length(all_MMs)
    name = all_MMs{kk};
    idx          = find(strcmp(basisSetIn.name,name));
    if ~isempty(idx)
        basisSetOut.name{actBasisFnct} = basisSetIn.name{idx};
        basisSetOut.fids(:,actBasisFnct,:)   = basisSetIn.fids(:,idx,:);
        basisSetOut.specs(:,actBasisFnct,:)  = basisSetIn.specs(:,idx,:);
        actBasisFnct = actBasisFnct + 1;
    end
end

basisSetOut.nMM = actBasisFnct - basisSetOut.nMets-1;

try
    basisSetOut.name(actBasisFnct:end) = [];
    basisSetOut.fids(:,actBasisFnct:end,:)   = [];
    basisSetOut.specs(:,actBasisFnct:end,:)  = [];
catch
end
basisSetOut.sz = size(basisSetOut.fids);

end