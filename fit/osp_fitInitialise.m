function [MRSCont] = osp_fitInitialise(MRSCont)
%% [MRSCont] = osp_fitInitialise(MRSCont)
%   This function initialises default basis sets and decides which
%   metabolites to include in the modeling process carried out by OspreyFit.
%
%   USAGE:
%       [MRSCont] = osp_fitInitialise(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-24)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-02-24: First version of the code.

% Find the right basis set (provided as *.mat file in Osprey basis set
% format)
if MRSCont.flags.isUnEdited
    % Extract TE from first dataset
    te = num2str(MRSCont.raw{1}.te);
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/philips/press' te '/basis_philips_press' te '.mat']); 
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/ge/press' te '/basis_ge_press' te '.mat']); 
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/siemens/press' te '/basis_siemens_press' te '.mat']); 
    end
elseif MRSCont.flags.isMEGA
    % Extract edit target from MRSCont
    editTarget = lower(MRSCont.opts.editTarget{1});
    % Extract TE from first dataset
    te = num2str(MRSCont.raw{1}.te);
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/philips/megapress_' editTarget te '/basis_philips_megapress_' editTarget te '.mat']);
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/ge/megapress_' editTarget te '/basis_ge_megapress_' editTarget te '.mat']);
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/siemens/megapress_' editTarget te '/basis_siemens_megapress_' editTarget te '.mat']);
    end
elseif MRSCont.flags.isHERMES
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hermes/basis_siemens_hermes.mat');
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/GE/HERMES/BASIS_noMM.mat');
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hermes/basis_siemens_hermes.mat');
    end
elseif MRSCont.flags.isHERCULES
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/philips/hercules-press/basis_philips_hercules-press.mat');
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/GE/HERCULES/BASIS.mat');
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hercules-press/basis_siemens_hercules-press.mat');
    end
end

% Clear existing basis set
MRSCont.fit.basisSet = [];
% Load the specified basis set
basisSet = load(MRSCont.opts.fit.basisSetFile);
basisSet = basisSet.BASIS;

% Generate the list of basis functions that are supposed to be included in
% the basis set
includeMetabs = MRSCont.opts.fit.includeMetabs;
metabList = fit_createMetabList(includeMetabs);
% Collect MMfit flag from the options determined in the job file
fitMM = MRSCont.opts.fit.fitMM;
% Create the modified basis set
basisSet = fit_selectMetabs(basisSet, metabList, fitMM);

% Determine the scaling factor between data and basis set for each dataset
for kk = 1:MRSCont.nDatasets
    MRSCont.fit.scale{kk} = max(real(MRSCont.processed.A{kk}.specs)) / max(max(max(real(basisSet.specs))));
end

% Save the modified basis set
MRSCont.fit.basisSet = basisSet;

end