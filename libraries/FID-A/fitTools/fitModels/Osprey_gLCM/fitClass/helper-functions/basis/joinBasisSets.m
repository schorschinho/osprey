function [basisSet] = joinBasisSets(basisSet, basisSim)
% Georg Oeltzschner, Johns Hopkins University 2023
% Joins two basis sets (intended to join a metabolite basis set with
% another basis set that contains simulated/parametrized MM/Lip signals).
% Untested for two arbitrary basis sets.

% Append basis struct
if size(basisSet.fids,3) == 1
    basisSet.fids   = horzcat(basisSet.fids, basisSim.fids);
    basisSet.specs  = horzcat(basisSet.specs, basisSim.specs);
    basisSet.nMets  = basisSet.nMets + basisSim.nMets;
    basisSet.nMM    = basisSet.nMM + basisSim.nMM;
    basisSet.sz     = size(basisSet.fids);
    basisSet.name   = horzcat(basisSet.name, basisSim.name);
else
    basisSet.fids   = horzcat(basisSet.fids, repmat(basisSim.fids,[1,1,basisSet.nExtra]));
    basisSet.specs  = horzcat(basisSet.specs, repmat(basisSim.specs,[1,1,basisSet.nExtra]));
    basisSet.nMets  = basisSet.nMets + basisSim.nMets;
    basisSet.nMM    = basisSet.nMM + basisSim.nMM;
    basisSet.sz     = size(basisSet.fids);
    basisSet.name   = horzcat(basisSet.name, basisSim.name);
end

end