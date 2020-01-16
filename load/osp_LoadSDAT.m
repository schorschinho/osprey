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

% Determine number of datasets
MRSCont.nDatasets = length(MRSCont.files);
if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets
        error('Number of specified reference files does not match number of specified metabolite files.');
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets
        error('Number of specified water files does not match number of specified metabolite files.');
    end
end

%% Get the data (loop over all datasets)
refLoadTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
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
fprintf('... done.\n');
toc(refLoadTime);

% Set flag
MRSCont.flags.coilsCombined     = 1;

end