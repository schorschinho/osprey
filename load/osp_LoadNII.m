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
        [~] = printLog('OspreyLoad', kk,1, MRSCont.nDatasets ,progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);   
        if ~(MRSCont.flags.didLoad == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

            % Read in the raw metabolite data.
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            raw 	= io_loadspec_niimrs(MRSCont.files{metab_ll,kk});
            % Add NIfTI-MRS information
            raw                           = osp_add_nii_mrs_field(raw,MRSCont.ver.Osp);
            % If the sequence flags are not parsed form the nii header
            if sum([raw.flags.isUnEdited raw.flags.isMEGA raw.flags.isHERMES raw.flags.isHERCULES]) == 0
                if MRSCont.flags.isUnEdited
                    raw.flags.isUnEdited = 1;
                end
                if MRSCont.flags.isMEGA
                    raw.flags.isMEGA = 1;
                end
                if MRSCont.flags.isHERMES
                    raw.flags.isHERMES = 1;
                end
                if MRSCont.flags.isHERCULES
                    raw.flags.isHERCULES = 1;
                end
            end

            % Read in the MM data.
            % Leave until we have example data.

            % Read in the raw water reference data.
            if MRSCont.flags.hasRef
                ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                raw_ref     = io_loadspec_niimrs(MRSCont.files_ref{ref_ll,kk});
                % If the sequence flags are not parsed form the nii header
                if sum([raw_ref.flags.isUnEdited raw_ref.flags.isMEGA raw_ref.flags.isHERMES raw_ref.flags.isHERCULES]) == 0
                    if MRSCont.flags.isUnEdited
                        raw_ref.flags.isUnEdited = 1;
                    end
                    if MRSCont.flags.isMEGA
                        raw_ref.flags.isMEGA = 1;
                    end
                    if MRSCont.flags.isHERMES
                        raw_ref.flags.isHERMES = 1;
                    end
                    if MRSCont.flags.isHERCULES
                        raw_ref.flags.isHERCULES = 1;
                    end
                end
                % Add NIfTI-MRS information
                raw_ref   = osp_add_nii_mrs_field(raw_ref,MRSCont.ver.Osp);
            end

            % Read in the short-TE water data.
            if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                raw_w       = io_loadspec_niimrs(MRSCont.files_w{w_ll,kk});
                % If the sequence flags are not parsed form the nii header
                if sum([raw_w.flags.isUnEdited raw_w.flags.isMEGA raw_w.flags.isHERMES raw_w.flags.isHERCULES]) == 0
                    if MRSCont.flags.isUnEdited
                        raw_w.flags.isUnEdited = 1;
                    end
                    if MRSCont.flags.isMEGA
                        raw_w.flags.isMEGA = 1;
                    end
                    if MRSCont.flags.isHERMES
                        raw_w.flags.isHERMES = 1;
                    end
                    if MRSCont.flags.isHERCULES
                        raw_w.flags.isHERCULES = 1;
                    end
                end
                % Add NIfTI-MRS information
                raw_w   = osp_add_nii_mrs_field(raw_w,MRSCont.ver.Osp);
            end

            
            


            % Set flag and save data under appropriate name
            if raw.dims.coils == 0
                MRSCont.flags.coilsCombined = 1;

                MRSCont.raw{metab_ll,kk}      = raw;

                if MRSCont.flags.hasRef
                    raw_ref = op_combine_water_subspecs(raw_ref,0);
                    MRSCont.raw_ref{ref_ll,kk}  = raw_ref;
                end

                if MRSCont.flags.hasWater
                    raw_w = op_combine_water_subspecs(raw_w,0);
                    MRSCont.raw_w{w_ll,kk}    = raw_w;
                end

            else
                MRSCont.flags.coilsCombined = 0;

                MRSCont.raw_uncomb{metab_ll,kk}      = raw;

                if MRSCont.flags.hasRef
                    MRSCont.raw_ref_uncomb{ref_ll,kk}  = raw_ref;
                end

                if MRSCont.flags.hasWater
                    MRSCont.raw_w_uncomb{w_ll,kk}    = raw_w;
                end

            end
        end
     end
    
    % Set flag and save data under appropriate name
    if raw.dims.coils == 0
        MRSCont.flags.coilsCombined = 1;
        % Do some houskeeping on the metabolite data if the nii-mrs conversion
        % did not work correctly 
        if (raw.flags.isMEGA || raw.flags.isHERMES || raw.flags.isHERCULES) && (raw.dims.subSpecs == 0)
            if raw.flags.isMEGA
                raw.fids = reshape(raw.fids,[raw.sz(raw.dims.t) raw.sz(raw.dims.averages)/2 2]);
                raw.specs = reshape(raw.specs,[raw.sz(raw.dims.t) raw.sz(raw.dims.averages)/2 2]);
                raw.averages = raw.sz(raw.dims.averages)/2;
                raw.sz = size(raw.fids);
                raw.subspecs = 2;
                raw.rawSubspecs = 2;
                raw.dims.subSpecs = 3;
            end
            if raw.flags.isHERMES || raw.flags.isHERCULES
                raw.fids = reshape(raw.fids,[raw.sz(raw.dims.t) raw.sz(raw.dims.averages)/4 4]);
                raw.specs = reshape(raw.specs,[raw.sz(raw.dims.t) raw.sz(raw.dims.averages)/4 4]);
                raw.averages = raw.sz(raw.dims.averages)/4;
                raw.sz = size(raw.fids);
                raw.subspecs = 4;
                raw.rawSubspecs = 4;  
                raw.dims.subSpecs = 3;
            end
        end
        MRSCont.raw{kk}      = raw;
        if MRSCont.flags.hasRef
            MRSCont.raw_ref{kk}  = raw_ref;
        end
        
        if MRSCont.flags.hasWater
            MRSCont.raw_w{kk}    = raw_w;
        end
    else
        MRSCont.flags.coilsCombined = 0;
        % Do some houskeeping on the metabolite data if the nii-mrs conversion
        % did not work correctly 
        if (raw.flags.isMEGA || raw.flags.isHERMES || raw.flags.isHERCULES) && (raw.dims.subSpecs == 0)
            if raw.flags.isMEGA
                raw.fids = reshape(raw.fids,[raw.sz(raw.dims.t) raw.sz(raw.dims.coils) raw.sz(raw.dims.averages)/2 2]);
                raw.specs = reshape(raw.specs,[raw.sz(raw.dims.t) raw.sz(raw.dims.coils) raw.sz(raw.dims.averages)/2 2]);
                raw.averages = raw.sz(raw.dims.averages)/2;
                raw.sz = size(raw.fids);
                raw.subspecs = 2;
                raw.rawSubspecs = 2;
                raw.dims.subSpecs = 4;
            end
            if raw.flags.isHERMES || raw.flags.isHERCULES
                raw.fids = reshape(raw.fids,[raw.sz(raw.dims.t) raw.sz(raw.dims.coils) raw.sz(raw.dims.averages)/4 4]);
                raw.specs = reshape(raw.specs,[raw.sz(raw.dims.t) raw.sz(raw.dims.coils) raw.sz(raw.dims.averages)/4 4]);
                raw.averages = raw.sz(raw.dims.averages)/4;
                raw.sz = size(raw.fids);
                raw.subspecs = 4;
                raw.rawSubspecs = 4;   
                raw.dims.subSpecs = 4;
            end
        end
        MRSCont.raw_uncomb{kk}      = raw;
        if MRSCont.flags.hasRef
            MRSCont.raw_ref_uncomb{kk}  = raw_ref;
        end
        
        if MRSCont.flags.hasWater
            MRSCont.raw_w_uncomb{kk}    = raw_w;
        end
    end

end
% Try to get some params from NIFTI
if isfield(raw.nii_mrs.hdr_ext,'Manufacturer')
    switch lower(raw.nii_mrs.hdr_ext.Manufacturer)
        case 'siemens'
                MRSCont.vendor = 'Siemens';
        case 'philips'
            MRSCont.vendor = 'Philips';
        case 'ge'
            MRSCont.vendor = 'GE';
    end
else
    MRSCont.vendor = 'NIFTI';
end
if isfield(raw.nii_mrs.hdr_ext,'ProtocolName') && ~isempty(raw.nii_mrs.hdr_ext.ProtocolName)
    MRSCont.seq = raw.nii_mrs.hdr_ext.ProtocolName;
else
    MRSCont.seq = 'NIFTI'; 
end

time = toc(refLoadTime);
[~] = printLog('done',time,MRSCont.nDatasets(1),MRSCont.nDatasets(2),progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 

MRSCont.runtime.Load = time;

end