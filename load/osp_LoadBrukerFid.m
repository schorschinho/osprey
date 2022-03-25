function [MRSCont] = osp_LoadBrukerFid(MRSCont)
%% [MRSCont] = osp_LoadBrukerFid(MRSCont)
%   This function reads raw data from Bruker FID datasets.
%
%   USAGE:
%       [MRSCont] = osp_LoadBrukerFid(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-07-06)
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

% Check that all datasets are complete
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
        if ~(MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

            % Read in the raw metabolite data.
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            raw = io_loadspec_bruk(MRSCont.files{metab_ll,kk});
            % Sequence flags
            if MRSCont.flags.isUnEdited
                raw.flags.isUnEdited = 1;
            end
            if MRSCont.flags.isMEGA
                raw.flags.isMEGA = 1;
            end
            MRSCont.raw{metab_ll,kk}      = raw;

            % Read in the raw MM data. re_mm
            if MRSCont.flags.hasMM %re_mm
                temp_ll = MRSCont.opts.MultipleSpectra.mm(ll);
                if MRSCont.flags.isUnEdited %re_mm
                    raw_mm = io_loadspec_bruk(MRSCont.files_mm{temp_ll,kk},1); %re_mm
                    [raw_mm] = op_rmempty(raw_mm); %re_mm
                    raw_mm.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA %re_mm
                    raw_mm = io_loadspec_bruk(MRSCont.files_mm{temp_ll,kk},2); %re_mm
                    if raw_mm.subspecs > 1 %re_mm
                        raw_mm_A               = op_takesubspec(raw_mm,1); %re_mm
                        [raw_mm_A]             = op_rmempty(raw_mm_A);            % Remove empty linesv
                        raw_mm_B               = op_takesubspec(raw_mm,2); %re_mm
                        [raw_mm_B]             = op_rmempty(raw_mm_B);            % Remove empty lines %re_mm
                        raw_mm                 = op_concatAverages(raw_mm_A,raw_mm_B); %re_mm
                    else %re_mm
                        [raw_mm] = op_rmempty(raw_mm); %re_mm
                    end %re_mm
                    raw_mm.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES%re_mm
                    raw_mm = io_loadspec_bruk(MRSCont.files_mm{temp_ll,kk},1); %re_mm
                    if raw_mm.subspecs > 1 %re_mm
                        raw_mm_A               = op_takesubspec(raw_mm,1); %re_mm
                        [raw_mm_A]             = op_rmempty(raw_mm_A);            % Remove empty lines %re_mm
                        raw_mm_B               = op_takesubspec(raw_mm,2); %re_mm
                        [raw_mm_B]             = op_rmempty(raw_mm_B);            % Remove empty lines %re_mm
                        raw_mm_C               = op_takesubspec(raw_mm,3); %re_mm
                        [raw_mm_C]             = op_rmempty(raw_mm_C);            % Remove empty lines %re_mm
                        raw_mm_D               = op_takesubspec(raw_mm,4); %re_mm
                        [raw_mm_D]             = op_rmempty(raw_mm_D);            % Remove empty lines %re_mm
                        raw_mm                 = op_concatAverages(raw_mm_A,raw_mm_B,raw_mm_C,raw_mm_D);   %re_mm
                    else %re_mm
                        [raw_mm] = op_rmempty(raw_mm);  %re_mm
                    end %re_mm
                    raw_mm.flags.isHERMES = 1;
                elseif MRSCont.flags.isHERCULES %re_mm
                    raw_mm = io_loadspec_bruk(MRSCont.files_mm{temp_ll,kk},1); %re_mm
                    if raw_mm.subspecs > 1 %re_mm
                        raw_mm_A               = op_takesubspec(raw_mm,1); %re_mm
                        [raw_mm_A]             = op_rmempty(raw_mm_A);            % Remove empty lines %re_mm
                        raw_mm_B               = op_takesubspec(raw_mm,2); %re_mm
                        [raw_mm_B]             = op_rmempty(raw_mm_B);            % Remove empty lines %re_mm
                        raw_mm_C               = op_takesubspec(raw_mm,3); %re_mm
                        [raw_mm_C]             = op_rmempty(raw_mm_C);            % Remove empty lines %re_mm
                        raw_mm_D               = op_takesubspec(raw_mm,4); %re_mm
                        [raw_mm_D]             = op_rmempty(raw_mm_D);            % Remove empty lines %re_mm
                        raw_mm                 = op_concatAverages(raw_mm_A,raw_mm_B,raw_mm_C,raw_mm_D);   %re_mm
                    else %re_mm
                        [raw_mm] = op_rmempty(raw_mm);  %re_mm
                    end %re_mm
                    raw_mm.flags.isHERCULES = 1;
                end
                MRSCont.raw_mm{temp_ll,kk}  = raw_mm;
            end %re_mm      

            % Read in the raw reference data.
            if MRSCont.flags.hasRef
                ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                if MRSCont.flags.isUnEdited
                    raw_ref = io_loadspec_bruk(MRSCont.files_ref{ref_ll,kk},1);
                    [raw_ref] = op_rmempty(raw_ref);
                    raw_ref.flags.isUnEdited = 1;                
                elseif MRSCont.flags.isMEGA
                    raw_ref = io_loadspec_bruk(MRSCont.files_ref{ref_ll,kk},2);
                    if raw_ref.subspecs > 1
                        raw_ref_A               = op_takesubspec(raw_ref,1);
                        [raw_ref_A]             = op_rmempty(raw_ref_A);            % Remove empty lines
                        raw_ref_B               = op_takesubspec(raw_ref,2);
                        [raw_ref_B]             = op_rmempty(raw_ref_B);            % Remove empty lines
                        raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B);
                    else
                        [raw_ref] = op_rmempty(raw_ref);
                    end
                    raw_ref.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES
                    raw_ref = io_loadspec_bruk(MRSCont.files_ref{ref_ll,kk},1);
                    if raw_ref.subspecs > 1
                        raw_ref_A               = op_takesubspec(raw_ref,1);
                        [raw_ref_A]             = op_rmempty(raw_ref_A);            % Remove empty lines
                        raw_ref_B               = op_takesubspec(raw_ref,2);
                        [raw_ref_B]             = op_rmempty(raw_ref_B);            % Remove empty lines
                        raw_ref_C               = op_takesubspec(raw_ref,3);
                        [raw_ref_C]             = op_rmempty(raw_ref_C);            % Remove empty lines
                        raw_ref_D               = op_takesubspec(raw_ref,4);
                        [raw_ref_D]             = op_rmempty(raw_ref_D);            % Remove empty lines
                        raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B,raw_ref_C,raw_ref_D);  
                    else
                        [raw_ref] = op_rmempty(raw_ref); 
                    end
                    raw_ref.flags.isHERMES = 1;
                elseif MRSCont.flags.isHERCULES
                    raw_ref = io_loadspec_bruk(MRSCont.files_ref{ref_ll,kk},1);
                    if raw_ref.subspecs > 1
                        raw_ref_A               = op_takesubspec(raw_ref,1);
                        [raw_ref_A]             = op_rmempty(raw_ref_A);            % Remove empty lines
                        raw_ref_B               = op_takesubspec(raw_ref,2);
                        [raw_ref_B]             = op_rmempty(raw_ref_B);            % Remove empty lines
                        raw_ref_C               = op_takesubspec(raw_ref,3);
                        [raw_ref_C]             = op_rmempty(raw_ref_C);            % Remove empty lines
                        raw_ref_D               = op_takesubspec(raw_ref,4);
                        [raw_ref_D]             = op_rmempty(raw_ref_D);            % Remove empty lines
                        raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B,raw_ref_C,raw_ref_D);  
                    else
                        [raw_ref] = op_rmempty(raw_ref); 
                    end
                    raw_ref.flags.isHERCULES = 1;
                end
                MRSCont.raw_ref{ref_ll,kk}  = raw_ref;
            end
            if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                raw_w   = io_loadspec_bruk(MRSCont.files_w{w_ll,kk},1);
                raw_w.flags.isUnEdited = 1;
                MRSCont.raw_w{w_ll,kk}    = raw_w;
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