function [MRSCont] = osp_LoadNII(MRSCont)
%% [MRSCont] = osp_LoadNII(MRSCont)
%   This function reads raw data from NIfTI-MRS files.
%
%   USAGE:
%       [MRSCont] = osp_LoadNII(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-07-21)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)


% Close any remaining open figures
close all;
warning('off','all');

%re_mm adding functionality to load MM data
if MRSCont.flags.hasMM 
    if ((length(MRSCont.files_mm) == 1) && (MRSCont.nDatasets>1))   %re_mm seems like specificy one MM file for a batch is also an option to plan to accomodate
        for kk=2:MRSCont.nDatasets %re_mm 
            MRSCont.files_mm{kk} = MRSCont.files_mm{1}; % re_mm allowable to specify one MM file for the whole batch
        end %re_mm 
    end   %re_mm 
    if ((length(MRSCont.files_mm) ~= MRSCont.nDatasets) )   %re_mm 
        msg = 'Number of specified MM files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end   %re_mm 
end   %re_mm 

if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets
        msg = 'Number of specified reference files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end
end

if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets
        msg = 'Number of specified water files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end
end

%% Get the data (loop over all datasets)
refLoadTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyLoad', kk, MRSCont.nDatasets ,progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);   
    if ~(MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

        % Read in the raw metabolite data.
        raw 	= io_loadspec_niimrs(MRSCont.files{kk});

        % Read in the MM data.
        % Leave until we have example data.
    
        % Read in the raw water reference data.
        if MRSCont.flags.hasRef
            raw_ref     = io_loadspec_niimrs(MRSCont.files_ref{kk});
        end
        
        % Read in the short-TE water data.
        if MRSCont.flags.hasWater
            raw_w       = io_loadspec_niimrs(MRSCont.files_w{kk});
        end
    end
end
time = toc(refLoadTime);
[~] = printLog('done', time, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI); 

% Set flag and save data under appropriate name
if raw.dims.coils == 0
    MRSCont.flags.coilsCombined = 1;
    
    MRSCont.raw{kk}      = raw;
    
    if MRSCont.flags.hasRef
        MRSCont.raw_ref{kk}  = raw_ref;
    end
    
    if MRSCont.flags.hasWater
        MRSCont.raw_w{kk}    = raw_w;
    end
    
else
    MRSCont.flags.coilsCombined = 0;
    
    MRSCont.raw_uncomb{kk}      = raw;
    
    if MRSCont.flags.hasRef
        MRSCont.raw_ref_uncomb{kk}  = raw_ref;
    end
    
    if MRSCont.flags.hasWater
        MRSCont.raw_w_uncomb{kk}    = raw_w;
    end
    
end

MRSCont.runtime.Load = time;

end