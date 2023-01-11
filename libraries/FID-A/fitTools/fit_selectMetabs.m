% fit_selectMetabs.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% basisSetOut = fit_selectMetabs(basisSetIn, metabList, fitMM)
% 
% DESCRIPTION:
% This function loads the basis set functions for the metabolites
% that have a positive flag set in the structure metabList, which has
% previously been created with fit_createMetabList.
%
% Only basis functions for the metabolites selected in fit_createMetabList
% will be included. This function will output a modified FID-A basis
% set container, which will be subsequently passed on to the fitting
% algorithm.
%
% The function will also include MM and lipid spectra, if present in
% the basis set, and if they were selected to be included.
%
% The other basis functions will be discarded.
% 
% OUTPUTS:
% basisSetOut = FID-A basis set container that only contains the basis
%               functions specified in metabList
%
% INPUTS:
% basisSetIn  = FID-A basis set container (created with fit_makeBasis).
% metabList   = Structure (created with fit_createMetabList) containing
%               flags for each metabolite/MM/lipid basis function that is 
%               supposed to be included in the basis set.
% fitMM       = Flag determining whether MM/lipid basis functions are being
%               kept in the output basis set. 1 = Yes (default), 0 = No.

function basisSetOut = fit_selectMetabs(basisSetIn, metabList, fitMM)

% Parse input arguments
if nargin<3
    fitMM = 1;
end
% Save all available metabolite names in a cell
all_mets = {'Ala','Asc','Asp','bHB','bHG','Cit','Cr','Cystat','CrCH2','EA','EtOH','fCho','GABA','GPC','GSH','Glc','Gln' ...
    ,'Glu','Gly','HCar','H2O','mI','Lac','NAA','NAAG','PCh','PCr','PE','Phenyl' ...
    ,'sI','Ser','Tau','Tyros','NAA_Ace','NAA_Asp'};

% Duplicate the input basis set
basisSetOut = basisSetIn;

% Check which metabolites are available in the basis set
metsInBasisSet = basisSetIn.name;
[metsToKeep,~,~] = intersect(metsInBasisSet, all_mets, 'stable');

% Check for each remaining metabolite in the basis set whether it has been 
% flagged to be included in the basis set in osp_FitInitialise.m:
idx_toKeep = zeros(length(metsToKeep),1);
for kk = 1:length(metsToKeep)
    if ~isfield(metabList, metsToKeep{kk}) || ~metabList.(metsToKeep{kk})
        idx_toKeep(kk) = 0;
    else
        idx_toKeep(kk) = 1;
    end
end
basisSetOut.nMets = sum(idx_toKeep);

% If the flag for including MM/lipid basis functions is set in osp_FitInitialise,
% include them now

all_MMs = {'MM09','MM12','MM14','MM17','MM20','MM22', 'MM27','MM30','MM32','Lip09','Lip13','Lip20','MM37','MM38','MM40','MM42','MMexp','MM_PRESS_PCC','MM_PRESS_CSO'};

% Check which of these are available in the basis set
MMsInBasisSet = basisSetIn.name;
[MMsToKeep,~,~] = intersect(MMsInBasisSet, all_MMs, 'stable');
idx_toKeepMM = zeros(length(MMsToKeep),1);
if fitMM
    for kk = 1:length(MMsToKeep)
        if ~isfield(metabList, MMsToKeep{kk}) || ~metabList.(MMsToKeep{kk})
            idx_toKeepMM(kk) = 0;
        else
            idx_toKeepMM(kk) = 1;
        end
    end
end
% if fitMM ==1
%     idx_toKeepMM = ones(length(MMsToKeep),1);
% else if fitMM == 2
%         idx_toKeepMM = zeros(8,1);
%         idx_toKeepMM(end) = 1;
%     else
%         idx_toKeepMM = zeros(length(MMsToKeep),1);    
%     end
% end

idx_toKeep = [idx_toKeep; idx_toKeepMM];
basisSetOut.nMM = sum(idx_toKeepMM);


% If the flag is set to zero, remove the name, the FIDs and the specs
% from the basis set. Also update the numbers for metabolite and MM/lipid
% basis functions
basisSetOut.name   = basisSetOut.name(logical(idx_toKeep));
basisSetOut.fids   = basisSetOut.fids(:,logical(idx_toKeep),:);
basisSetOut.specs  = basisSetOut.specs(:,logical(idx_toKeep),:);
basisSetOut.sz = size(basisSetOut.fids);


end