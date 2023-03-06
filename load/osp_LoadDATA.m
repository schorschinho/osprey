function [MRSCont] = osp_LoadDATA(MRSCont)
%% [MRSCont] = osp_LoadSDAT(MRSCont)
%   Reads Philips DATA/LIST files it is adapted from Gannet.
%
%   Author:
%       Dr.Helge Zoellner (Johns Hopkins University, 2020-10-02)
%       hzoelln2@jhmi.edu
%   
%   Credits:
%       This code uses the function
%       loadRawKspace.m
%       from the excellent "Matlab raw kspace tools" toolbox
%       (Wouter Potters, Academic Medical Center, Amsterdam, NL)
%       https://bitbucket.org/wpotters/matlab-raw-kspace-tools
%
%   History:
%       2020-10-02: First version.


% Close any remaining open figures
close all;
warning('off','all');
if MRSCont.flags.hasMM %re_mm adding functionality to load MM data
    if ((length(MRSCont.files_mm) == 1) && (MRSCont.nDatasets(1)>1))   %re_mm seems like specificy one MM file for a batch is also an option to plan to accomodate
        for kk=2:MRSCont.nDatasets(1) %re_mm 
            MRSCont.files_mm{kk} = MRSCont.files_mm{1}; % re_mm allowable to specify one MM file for the whole batch
        end %re_mm 
    end   %re_mm 
    if ((length(MRSCont.files_mm) ~= MRSCont.nDatasets(1)) )   %re_mm 
        msg = 'Number of specified MM files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end   %re_mm 
end   %re_mm 
if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets(1)
        msg = 'Number of specified reference files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets(1)
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
for kk = 1:MRSCont.nDatasets(1)
    for ll = 1: 1:MRSCont.nDatasets(2)
        [~] = printLog('OspreyLoad',kk,ll,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);  

        if ~(MRSCont.flags.didLoad == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

            % Read in the raw metabolite data. Since the Philips DATA loader needs
            % to know the number of sub-spectra (e.g. from spectral editing), the
            % type of sequence needs to be differentiated here already.
            if MRSCont.flags.hasStatfile
                statFile = MRSCont.file_stat;
            else
                statFile = [];
            end
            
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            if MRSCont.flags.isUnEdited
                [raw,raw_ref]=io_loadspec_data(MRSCont.files{metab_ll,kk},1,kk,statFile);
                raw.flags.UnEdited = 1;
                raw_ref.flags.UnEdited = 1;
            elseif MRSCont.flags.isMEGA
                [raw,raw_ref]=io_loadspec_data(MRSCont.files{metab_ll,kk},2,kk,statFile);
                raw.flags.isMEGA = 1;
                raw_ref.flags.isMEGA = 1;
            elseif MRSCont.flags.isHERMES
                [raw,raw_ref]=io_loadspec_data(MRSCont.files{metab_ll,kk},4,kk,statFile);
                raw.flags.isHERMES = 1;
                raw_ref.flags.isHERMES = 1;
            elseif MRSCont.flags.isHERCULES
                [raw,raw_ref]=io_loadspec_data(MRSCont.files{metab_ll,kk},4,kk,statFile);
                raw.flags.isHERCULES = 1;
                raw_ref.flags.isHERCULES = 1;
            end
            % Add NIfTI-MRS information
            raw  = osp_add_nii_mrs_field(raw,MRSCont.ver.Osp);
            MRSCont.raw_uncomb{metab_ll,kk}      = raw;

            if ~MRSCont.flags.hasRef && ~isempty(raw_ref)
                % Add NIfTI-MRS information
                raw_ref  = osp_add_nii_mrs_field(raw_ref,MRSCont.ver.Osp);
                MRSCont.raw_ref_uncomb{metab_ll,kk}  = raw_ref;
                if kk == MRSCont.nDatasets
                    MRSCont.flags.hasRef = 1;
                    MRSCont.opts.MultipleSpectra.ref = MRSCont.opts.MultipleSpectra.metab;
                end
            else if MRSCont.flags.hasRef
                    ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                    [~,raw_ref]=io_loadspec_data(MRSCont.files_ref{ref_ll,kk},1,kk,statFile);
                    % Add NIfTI-MRS information
                    raw_ref  = osp_add_nii_mrs_field(raw_ref,MRSCont.ver.Osp);
                    MRSCont.raw_ref_uncomb{ref_ll,kk}  = raw_ref;
                end
            end

            if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                [~,raw_w]=io_loadspec_data(MRSCont.files_w{w_ll,kk},1,kk,statFile);
                raw_w.flags.UnEdited = 1;
                % Add NIfTI-MRS information
                raw_w  = osp_add_nii_mrs_field(raw_w,MRSCont.ver.Osp);
                MRSCont.raw_w_uncomb{w_ll,kk}    = raw_w;
            end
        end
    end
    
end


time = toc(refLoadTime);
[~] = printLog('done',time,MRSCont.nDatasets(1),MRSCont.nDatasets(2),progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
% Set flag
MRSCont.flags.coilsCombined     = 0;
MRSCont.runtime.Load = time;
end
