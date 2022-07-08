function [MRSCont] = osp_LoadDICOM(MRSCont)
%% [MRSCont] = osp_LoadDICOM(MRSCont)
%   This function reads raw (but coil-combined) data from DICOM files
%   produced by Siemens sequences (*.dcm or *.ima).
%
%   USAGE:
%       [MRSCont] = osp_LoadDICOM(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-19)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-19: First version of the code.

% Close any remaining open figures
close all;
warning('off','all');
if MRSCont.flags.hasMM
    if length(MRSCont.files_mm) ~= MRSCont.nDatasets(1)
        msg = 'Number of specified metabolite-nulled files does not match number of specified metabolite files.';
        fprintf(msg);
        error(msg);
    end
end
if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets(1)
        msg = 'Number of specified reference files does not match number of specified metabolite files.';
        fprintf(msg);
        error(msg);
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets(1)
        msg = 'Number of specified water files does not match number of specified metabolite files.';
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
            % Read in the raw metabolite data.
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            raw = io_loadspec_dicom(MRSCont.files{metab_ll,kk});
            % Sequence flags
            if MRSCont.flags.isUnEdited
                raw.flags.isUnEdited = 1;
            end
            if MRSCont.flags.isMEGA
                raw.flags.isMEGA = 1;
            end
            MRSCont.raw{ll,kk}      = raw;

            % Read in the raw reference data.
            if MRSCont.flags.hasRef
                ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                raw_ref = io_loadspec_dicom(MRSCont.files_ref{ref_ll,kk});
                raw_ref = op_combine_water_subspecs(raw_ref,0);
                MRSCont.raw_ref{ref_ll,kk}  = raw_ref;
                end
                if MRSCont.flags.isUnEdited
                    MRSCont.raw_ref{ref_ll,kk}.flags.isUnEdited = 1;
                end
                if MRSCont.flags.isMEGA
                    MRSCont.raw_ref{ref_ll,kk}.flags.isMEGA = 1;
                end
            end
            if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                raw_w   = io_loadspec_dicom(MRSCont.files_w{w_ll,kk});
                raw_w = op_combine_water_subspecs(raw_w,0);
                MRSCont.raw_w{w_ll,kk}    = raw_w;
                if MRSCont.flags.isUnEdited
                    MRSCont.raw_w{w_ll,kk}.flags.isUnEdited = 1;
                end
                if MRSCont.flags.isMEGA
                    MRSCont.raw_w{w_ll,kk}.flags.isMEGA = 1;
                end
            end
            if MRSCont.flags.hasMM
                temp_ll = MRSCont.opts.MultipleSpectra.mm(ll);
                raw_mm   = io_loadspec_dicom(MRSCont.files_mm{temp_ll,kk});
                MRSCont.raw_mm{temp_ll,kk}    = raw_mm;
                if MRSCont.flags.isUnEdited
                    MRSCont.raw_mm{temp_ll,kk}.flags.isUnEdited = 1;
                end
                if MRSCont.flags.isMEGA
                    MRSCont.raw_mm{temp_ll,kk}.flags.isMEGA = 1;
                end
            end        
        end
    end
end
time = toc(refLoadTime);
[~] = printLog('done',time,ll,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
% Set flag
MRSCont.flags.coilsCombined     = 1;
MRSCont.runtime.Load = time;
end