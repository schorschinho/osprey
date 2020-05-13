function [MRSCont] = osp_fitWater(MRSCont, kk, fitWhich)
%% [MRSCont] = LCG_fitWater(MRSCont)
%   This function initializes and runs fitting of the water reference and
%   short-TE water scans.
%
%   USAGE:
%       [MRSCont] = osp_fitWater(MRSCont, fitWhich);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = Index of the cell containing the spectrum to be fit.
%       fitWhich    = String determining whether this function is performed
%                     on the water reference or the short-TE water scan.
%                     OPTIONS: - 'ref', 'w'
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-09)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-09: First version of the code.


% Find the basis set index corresponding to water
basisSet            = MRSCont.fit.basisSet;
h2o_idx = find(strcmp(basisSet.name, 'H2O'));
idx_toKeep = zeros(basisSet.nMets + basisSet.nMM,1);
idx_toKeep(h2o_idx) = 1;

% Remove the name, the FIDs and the specs of everything but water
% from the basis set.
basisSet.name   = basisSet.name(logical(idx_toKeep));
basisSet.fids   = basisSet.fids(:,logical(idx_toKeep),1); % index 1 because this is a GSH-OFF spectrum in edited data
basisSet.specs  = basisSet.specs(:,logical(idx_toKeep),1);
basisSet.nMets  = 1; basisSet.nMM = 0;

%%  Construct the basis functions and the spectrum that is to be fit.
% Apply scaling factor to the data
dataToFit       = MRSCont.processed.(fitWhich){kk};
dataToFit       = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
% Extract fit options
fitOpts     = MRSCont.opts.fit;
fitModel    = fitOpts.method;

% Call the fit function
[fitParamsWater, resBasisSetWater]  = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);
        
%% Save back the fit parameters to MRSCont
% Discern whether the reference or the short-TE water scan was selected
switch fitWhich
    case 'ref'
        str = 'ref';
    case 'w'
        str = 'w';
end
% Write
MRSCont.fit.resBasisSet.(str).water{kk}      = resBasisSetWater;
MRSCont.fit.results.(str).fitParams{kk}   = fitParamsWater;
        
end