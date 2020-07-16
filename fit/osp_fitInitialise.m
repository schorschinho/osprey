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

ext = 0;

if strcmp(MRSCont.vendor,'GE') % Still need to find a way to destinguish GE sequences
    seq = 'press';
end
if contains(seq,'gaba_par')
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
    editTarget = lower(MRSCont.opts.editTarget{1});  
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/philips/mega/' seq '/' editTarget te '/basis_philips_megapress_' editTarget te '.mat']);
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/ge/mega/' seq '/' editTarget te '/basis_ge_megapress_' editTarget te '.mat']);
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which(['fit/basissets/siemens/mega/' seq '/' editTarget te '/basis_siemens_megapress_' editTarget te '.mat']);
    end
elseif MRSCont.flags.isHERMES
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hermes/basis_siemens_hermes.mat');
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hermes/basis_siemens_hermes.mat');
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hermes/basis_siemens_hermes.mat');
    end
elseif MRSCont.flags.isHERCULES
    switch MRSCont.vendor
        case 'Philips'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/philips/hercules-press/basis_philips_hercules-press.mat');
        case 'GE'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/philips/hercules-press/basis_philips_hercules-press.mat');
        case 'Siemens'
            MRSCont.opts.fit.basisSetFile        = which('fit/basissets/siemens/hercules-press/basis_siemens_hercules-press.mat');
    end
end

% Clear existing basis set
MRSCont.fit.basisSet = [];

% Check if automated basis set pick worked, otherwise the basis set from
% the user folder is loaded.
if isempty(MRSCont.opts.fit.basisSetFile)
    addpath( which('fit/basissets'));
    MRSCont.opts.fit.basisSetFile = which('fit/basissets/user/BASIS_MM.mat');
    if isempty(MRSCont.opts.fit.basisSetFile)
        error('There is no appropriate basis set to model your data. Please supply a sufficient basis set in Osprey .mat format in the fit/basissets/user/BASIS_MM.mat file! ');
    else
        ext = 1;
    end
end

% Load the specified basis set or the user basis set file
basisSet = load(MRSCont.opts.fit.basisSetFile);
basisSet = basisSet.BASIS;

idx_H2O          = find(strcmp(basisSet.name,'H2O'));
if isempty(idx_H2O)
    WaterbasisSetFile        = which(['fit/basissets/philips/unedited/' seq '/' te '/basis_philips_' seq te '.mat']); 
    if isempty(WaterbasisSetFile)
        WaterbasisSetFile        = which(['fit/basissets/philips/unedited/press/' te '/basis_philips_press' te '.mat']);
    end
    if isempty(WaterbasisSetFile)
        WaterbasisSetFile        = which(['fit/basissets/philips/unedited/press/30/basis_philips_press30.mat']);
    end
    if isempty(MRSCont.opts.fit.basisSetFile)
        error('There is no appropriate water basis set to model your data. Please supply a sufficient basis set in Osprey .mat format in the fit/basissets/user/BASIS_MM.mat file! ');
    end
    basisSetWater = load(WaterbasisSetFile);
    basisSetWater = basisSetWater.BASIS;
    basisSet.name{basisSet.nMets+basisSet.nMM+1} = 'H2O';
    idx_H2O = find(strcmp(basisSetWater.name,'H2O'));
    basisSet.fids(:,basisSet.nMets+basisSet.nMM+1) = basisSetWater.fids(:,idx_H2O) .* basisSetWater.scale ./ basisSet.scale;
    basisSet.specs(:,basisSet.nMets+basisSet.nMM+1) = basisSetWater.specs(:,idx_H2O) .* basisSetWater.scale ./ basisSet.scale;
    basisSet.sz = size(basisSet.fids);
    basisSet.nMets = basisSet.nMets + 1;
    [~,maxH2Oind] = max(real(basisSet.specs(:,basisSet.nMets+basisSet.nMM)));
    ppmmaxH2O = basisSet.ppm(maxH2Oind);
    frqshift=(ppmmaxH2O-4.65)*basisSet.Bo*42.577;
    basisSet.fids(:,basisSet.nMets+basisSet.nMM)=basisSet.fids(:,basisSet.nMets+basisSet.nMM).*exp(-1i*basisSet.t*frqshift*2*pi)';
    basisSet.specs(:,basisSet.nMets+basisSet.nMM)=fftshift(fft(basisSet.fids(:,basisSet.nMets+basisSet.nMM),[],1),1);
end
    
% Generate the list of basis functions that are supposed to be included in
% the basis set
if ext
    % Sort basis set file according to Osprey conventions
    basisSet = fit_sortBasisSet(basisSet);
else  
    % To do: Interface with interactive user input
    metabList = fit_createMetabList(MRSCont.opts.fit.includeMetabs);
    % Collect MMfit flag from the options determined in the job file
    fitMM = MRSCont.opts.fit.fitMM;
    % Create the modified basis set
    basisSet = fit_selectMetabs(basisSet, metabList, fitMM);
end

% Determine the scaling factor between data and basis set for each dataset
for kk = 1:MRSCont.nDatasets
    MRSCont.fit.scale{kk} = max(real(MRSCont.processed.A{kk}.specs)) / max(max(max(real(basisSet.specs))));
end


% Save the modified basis set
MRSCont.fit.basisSet = basisSet;

end