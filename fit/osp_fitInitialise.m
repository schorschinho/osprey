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
% Extract TE and sequence from first dataset
te = num2str(MRSCont.raw{1}.te);
seq = lower(MRSCont.raw{1}.seq);
seq = seq(~ismember(seq, char([10 13]))); % remove return or carriage return

if strcmp(MRSCont.vendor,'GE') % Still need to find a way to destinguish GE sequences
    seq = 'press';
end
if contains(seq,'press')
    seq = 'press';
end
if contains(seq,'slaser')
    seq = 'slaser';
end
    
if MRSCont.flags.isUnEdited
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/philips/unedited/' seq '/' te '/basis_philips_' seq te '.mat']); 
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/ge/unedited/' seq '/' te '/basis_ge_' seq te '.mat']); 
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/siemens/unedited/' seq '/' te '/basis_siemens_' seq te '.mat']); 
    end
elseif MRSCont.flags.isMEGA
    switch MRSCont.vendor
        case 'Philips'
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/ge/mega/' seq '_' editTarget te '/basis_ge_megapress_' editTarget te '.mat']);
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/siemens/mega/' seq '_' editTarget te '/basis_siemens_megapress_' editTarget te '.mat']);
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
% To do: Interface with interactive user input
metabList = fit_createMetabList;
% Collect MMfit flag from the options determined in the job file
fitMM = MRSCont.opts.fit.fitMM;
% Create the modified basis set
basisSet = fit_selectMetabs(basisSet, metabList, fitMM);

% % Check that basis set and data have the same frequency axis orientation,
% % and flip it if necessary. For real-life MRS data, Osprey data structures
% % should always have the frequency axis running from high ppm to low ppm
% % values.
% polarityPPMData  = issorted(MRSCont.processed.A{1}.ppm, 'ascend');
% polarityPPMBasis = issorted(basisSet.ppm, 'ascend');
% if polarityPPMData == 1 && polarityPPMBasis == 0
%     % invert ppm axis
%     basisSet.ppm = wrev(basisSet.ppm);
%     % re-scale and fft FIDs
%     basisSet.specs = fftshift(ifft(conj(basisSet.fids),[],1),1) .* length(basisSet.ppm);
% end

% Determine the scaling factor between data and basis set for each dataset
for kk = 1:MRSCont.nDatasets
    MRSCont.fit.scale{kk} = max(real(MRSCont.processed.A{kk}.specs)) / max(max(max(real(basisSet.specs))));
end

% Save the modified basis set
MRSCont.fit.basisSet = basisSet;

end