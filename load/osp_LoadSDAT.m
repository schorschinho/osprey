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
    [~] = printLog('OspreyLoad',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);   
    if ~(MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

        % Read in the raw metabolite data. Since the Philips SDAT loader needs
        % to know the number of sub-spectra (e.g. from spectral editing), the
        % type of sequence needs to be differentiated here already.
        if MRSCont.flags.isUnEdited
            raw         = io_loadspec_sdat(MRSCont.files{kk},1);
        elseif MRSCont.flags.isMEGA
            raw         = io_loadspec_sdat(MRSCont.files{kk},2);
        elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
            raw         = io_loadspec_sdat(MRSCont.files{kk},4);
        end
        MRSCont.raw{kk}      = raw;
        
        % Read in the raw MM data. re_mm
        if MRSCont.flags.hasMM %re_mm
            if MRSCont.flags.isUnEdited %re_mm
                raw_mm = io_loadspec_sdat(MRSCont.files_mm{kk},1); %re_mm
                [raw_mm] = op_rmempty(raw_mm); %re_mm
            elseif MRSCont.flags.isMEGA %re_mm
                raw_mm = io_loadspec_sdat(MRSCont.files_mm{kk},2); %re_mm
                if raw_mm.subspecs > 1 %re_mm
                    raw_mm_A               = op_takesubspec(raw_mm,1); %re_mm
                    [raw_mm_A]             = op_rmempty(raw_mm_A);            % Remove empty linesv
                    raw_mm_B               = op_takesubspec(raw_mm,2); %re_mm
                    [raw_mm_B]             = op_rmempty(raw_mm_B);            % Remove empty lines %re_mm
                    raw_mm                 = op_concatAverages(raw_mm_A,raw_mm_B); %re_mm
                else %re_mm
                    [raw_mm] = op_rmempty(raw_mm); %re_mm
                end %re_mm
            elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES %re_mm
                raw_mm = io_loadspec_sdat(MRSCont.files_mm{kk},1); %re_mm
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
            end
            MRSCont.raw_mm{kk}  = raw_mm;
        end %re_mm      
    
        % Read in the raw reference data.
        if MRSCont.flags.hasRef
            if MRSCont.flags.isUnEdited
                raw_ref = io_loadspec_sdat(MRSCont.files_ref{kk},1);
                [raw_ref] = op_rmempty(raw_ref);
            elseif MRSCont.flags.isMEGA
                raw_ref = io_loadspec_sdat(MRSCont.files_ref{kk},2);
                if raw_ref.subspecs > 1
                    raw_ref_A               = op_takesubspec(raw_ref,1);
                    [raw_ref_A]             = op_rmempty(raw_ref_A);            % Remove empty lines
                    raw_ref_B               = op_takesubspec(raw_ref,2);
                    [raw_ref_B]             = op_rmempty(raw_ref_B);            % Remove empty lines
                    raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    [raw_ref] = op_rmempty(raw_ref);
                end
            elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                raw_ref = io_loadspec_sdat(MRSCont.files_ref{kk},1);
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
            end
            MRSCont.raw_ref{kk}  = raw_ref;
        end
        if MRSCont.flags.hasWater
            raw_w   = io_loadspec_sdat(MRSCont.files_w{kk},1);
            MRSCont.raw_w{kk}    = raw_w;
        end
    end
end
time = toc(refLoadTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
% Set flag
MRSCont.flags.coilsCombined     = 1;
MRSCont.runtime.Load = time;
end