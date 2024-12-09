function [MRSCont] = osp_LoadSDAT(MRSCont)
%% [MRSCont] = osp_LoadSDAT(MRSCont)
%   This function reads raw (but coil-combined) data from Philips SDAT files.
%
%   USAGE:
%       [MRSCont] = osp_LoadSDAT(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-20)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-20: First version of the code.

% Close any remaining open figures
close all;
warning('off','all');
if MRSCont.flags.hasMM %re_mm adding functionality to load MM data
    if ((length(MRSCont.files_mm) == 1) && (MRSCont.nDatasets(1)>1))   %re_mm seems like specificy one MM file for a batch is also an option to plan to accomodate
        for kk=2:MRSCont.nDatasets(1) %re_mm 
            MRSCont.files_mm{kk} = MRSCont.files_mm{1}; % re_mm allowable to specify one MM file for the whole batch
        end %re_mm 
    end   %re_mm 
    if ((size(MRSCont.files_mm,2) ~= MRSCont.nDatasets(1)) )   %re_mm 
        msg = 'Number of specified MM files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end   %re_mm 
end   %re_mm 
if MRSCont.flags.hasRef
    if size(MRSCont.files_ref,2) ~= MRSCont.nDatasets(1)
        msg = 'Number of specified reference files does not match number of specified metabolite files.'; %re_mm 
        fprintf(msg);
        error(msg);
    end
end
if MRSCont.flags.hasWater
    if size(MRSCont.files_w,2) ~= MRSCont.nDatasets(1)
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

            % Read in the raw metabolite data. Since the Philips SDAT loader needs
            % to know the number of sub-spectra (e.g. from spectral editing), the
            % type of sequence needs to be differentiated here already.
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            if MRSCont.flags.isUnEdited
                raw         = io_loadspec_sdat(MRSCont.files{metab_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                raw.flags.isUnEdited = 1;
            elseif MRSCont.flags.isMEGA
                raw         = io_loadspec_sdat(MRSCont.files{metab_ll,kk},2,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                raw.flags.isMEGA = 1;
            elseif MRSCont.flags.isHERMES
                raw         = io_loadspec_sdat(MRSCont.files{metab_ll,kk},4,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                raw.flags.isHERMES = 1;
            elseif MRSCont.flags.isHERCULES
                raw         = io_loadspec_sdat(MRSCont.files{metab_ll,kk},4,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                raw.flags.isHERCULES = 1;
            end
            % Add NIfTI-MRS information
            raw                           = osp_add_nii_mrs_field(raw,MRSCont.ver.Osp);
            MRSCont.raw{metab_ll,kk}      = raw;

            % Read in the raw MM data. re_mm
            if MRSCont.flags.hasMM %re_mm
                temp_ll = MRSCont.opts.MultipleSpectra.mm(ll);
                if MRSCont.flags.isUnEdited %re_mm                    
                    raw_mm = io_loadspec_sdat(MRSCont.files_mm{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle); %re_mm
                    [raw_mm] = op_rmempty(raw_mm); %re_mm
                    raw_mm.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA %re_mm
                    raw_mm = io_loadspec_sdat(MRSCont.files_mm{temp_ll,kk},2,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle); %re_mm
                    raw_mm.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES%re_mm
                    raw_mm = io_loadspec_sdat(MRSCont.files_mm{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle); %re_mm
                    raw_mm.flags.isHERMES = 1;
                elseif MRSCont.flags.isHERCULES %re_mm
                    raw_mm = io_loadspec_sdat(MRSCont.files_mm{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle); %re_mm
                    raw_mm.flags.isHERCULES = 1;
                end
                % Add NIfTI-MRS information
                raw_mm                      = osp_add_nii_mrs_field(raw_mm,MRSCont.ver.Osp);
                MRSCont.raw_mm{temp_ll,kk}  = raw_mm;
            end %re_mm      

            % Read in the raw reference data.
            if MRSCont.flags.hasRef
                ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                if MRSCont.flags.isUnEdited                    
                    raw_ref = io_loadspec_sdat(MRSCont.files_ref{ref_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_ref = op_combine_water_subspecs(raw_ref,0);
                    raw_ref.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA
                    try % GO 03042024: this works for the JHU patch only but throws an error for product MEGA
                        raw_ref = io_loadspec_sdat(MRSCont.files_ref{ref_ll,kk},2,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    catch % GO 03042024: this works for product MEGA
                        raw_ref = io_loadspec_sdat(MRSCont.files_ref{ref_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    end
                    raw_ref = op_combine_water_subspecs(raw_ref,0);
                    raw_ref.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES
                    raw_ref = io_loadspec_sdat(MRSCont.files_ref{ref_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_ref = op_combine_water_subspecs(raw_ref,0);
                    raw_ref.flags.isHERMES = 1;
                 elseif MRSCont.flags.isHERCULES
                    raw_ref = io_loadspec_sdat(MRSCont.files_ref{ref_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_ref = op_combine_water_subspecs(raw_ref,0);
                    raw_ref.flags.isHERCULES = 1;
                end
                % Add NIfTI-MRS information
                raw_ref                     = osp_add_nii_mrs_field(raw_ref,MRSCont.ver.Osp);
                MRSCont.raw_ref{ref_ll,kk}  = raw_ref;
            end
            if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                raw_w   = io_loadspec_sdat(MRSCont.files_w{w_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                raw_w = op_combine_water_subspecs(raw_w,0);
                raw_w.flags.isUnEdited = 1;
                % Add NIfTI-MRS information
                raw_w                     = osp_add_nii_mrs_field(raw_w,MRSCont.ver.Osp);
                MRSCont.raw_w{w_ll,kk}    = raw_w;
            end

            % Read in the MM reference data.
            if MRSCont.flags.hasMMRef
                temp_ll = MRSCont.opts.MultipleSpectra.mm_ref(ll);
                if MRSCont.flags.isUnEdited
                    raw_mm_ref = io_loadspec_sdat(MRSCont.files_mm_ref{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_mm_ref = op_combine_water_subspecs(raw_mm_ref,0);
                    raw_mm_ref.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA
                    raw_mm_ref = io_loadspec_sdat(MRSCont.files_mm_ref{temp_ll,kk},2,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_mm_ref = op_combine_water_subspecs(raw_mm_ref,0);
                    raw_mm_ref.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES
                    raw_mm_ref = io_loadspec_sdat(MRSCont.files_mm_ref{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_mm_ref = op_combine_water_subspecs(raw_mm_ref,0);
                    raw_mm_ref.flags.isHERMES = 1;
                 elseif MRSCont.flags.isHERCULES
                    raw_mm_ref = io_loadspec_sdat(MRSCont.files_mm_ref{temp_ll,kk},1,MRSCont.flags.isSERIES,MRSCont.opts.load.undoPhaseCycle);
                    raw_mm_ref = op_combine_water_subspecs(raw_mm_ref,0);
                    raw_mm_ref.flags.isHERCULES = 1;
                end
                % Add NIfTI-MRS information
                raw_mm_ref                      = osp_add_nii_mrs_field(raw_mm_ref,MRSCont.ver.Osp);
                MRSCont.raw_mm_ref{temp_ll,kk}  = raw_mm_ref;
            end
        end
    end
end
time = toc(refLoadTime);
[~] = printLog('done',time,MRSCont.nDatasets(1),MRSCont.nDatasets(2),progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
% Set flag
MRSCont.flags.coilsCombined     = 1;
MRSCont.runtime.Load = time;
end