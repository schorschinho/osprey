function [MRSCont] = osp_fitUnEdited(MRSCont)
%% [MRSCont] = osp_fitUnEdited(MRSCont)
%   This function performs spectral fitting of unedited MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitUnEdited(MRSCont, basisSet);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       basisSet    = Osprey basis set container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-12)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-12: First version of the code.


% Loop over all the datasets here
metFitTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('\nFitting metabolite spectra from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    % Apply scaling factor to the data
    dataToFit   = MRSCont.processed.A{kk};
    dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
    % Extract fit options
    fitOpts     = MRSCont.opts.fit;
    fitModel    = fitOpts.method;
    
    % Call the fit function
    basisSet            = MRSCont.fit.basisSet;
    [fitParams, resBasisSet] = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);
    
    % Save back the basis set and fit parameters to MRSCont
    MRSCont.fit.basisSet                    = basisSet;
    MRSCont.fit.resBasisSet.off{kk}             = resBasisSet;
    MRSCont.fit.results.off.fitParams{kk}   = fitParams;
    
    %% end time counter
    if isequal(kk, MRSCont.nDatasets)
        fprintf('... done.\n');
        toc(metFitTime);
    end
end


end