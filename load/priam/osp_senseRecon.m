function MRSCont = osp_senseRecon(MRSCont)
%% MRS_struct = osp_senseRecon(MRS_struct)
%   Reads Siemens TWIX files (*.dat) and removes participant information. New
%   de-identified TWIX files are then output, with filenames appended with
%   '_noID'. The original files are not overwritten.
%
%   Input:
%       TWIXDeIdentify, by itself, de-identifies all TWIX files found
%       within the current directory.
%
%   Output:
%       c = {'MRS_01.dat', 'MRS_02.dat'};
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-15)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-15: First version of the code.
%

%% Setup paths and filenames, and load coil reference scan
% Generate output folder for the Reconstruction
% Save the output structure to the output folder
outputFolder    = MRSCont.outputFolder;

if ~exist([outputFolder filesep 'SenseReconstruction'],'dir')
    mkdir([outputFolder filesep 'SenseReconstruction']);
end


%%
% Load coil reference scan
% If existing coil sensitivity matrices are saved in the SenseReconstruction
% folder, load them here. Otherwise, start loading them.

for kk = 1:MRSCont.nDatasets
    MRSCont.SENSE{kk}.sens_from_ref_scan = 0; % 1 not implemented yet, GO 11/01/2016
    % Ask for spatial separation of the voxels that was entered into the exam
    % card (in mm, positive value)
    MRSCont.SENSE{kk}.vox_sep = -60;
    % Get the input file name
        [path,filename,~]   = fileparts(MRSCont.files{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_split          = regexp(path,filesep,'split');
        if length(path_split) > 2
            saveName = [path_split{end-1} '_' path_split{end} '_' filename];
        end
        % Generate file name for the voxel reconstruction file to be saved under
        senseCoilFile            = fullfile(outputFolder ,'SenseReconstruction', [saveName '_ref_scan_sense_coil_img.mat']);
        volumeCoilFile            = fullfile(outputFolder ,'SenseReconstruction', [saveName '_ref_scan_volume_coil_img.mat']);
        MRSCont.SENSE{kk}.senseCoilFile = senseCoilFile;
        MRSCont.SENSE{kk}.volumeCoilFile = volumeCoilFile;
    if (exist(senseCoilFile,'file') ||  exist(volumeCoilFile,'file')) %fprintf('Found existing coil reference data. Loading...\n%s\n',[spec_path filesep 'GannetRecon_output' filesep 'ref_scan_sense_coil_img.mat']);
        load(senseCoilFile);
        load(volumeCoilFile);
        disp('Loading coil reference data finished!');
    else
        [cpx_file,Ref,img_array,img_body,noise_array,noise_body] = loadRefScan(MRSCont.files{kk});
        save(senseCoilFile,'cpx_file','Ref','img_array','noise_array','noise_body');
        save(volumeCoilFile,'cpx_file','Ref','img_body','noise_array','noise_body');
    end

    %% Calculate the actual unfolding matrix
    MRSCont = calcUnfoldingMatrix(MRSCont, cpx_file, Ref, img_array, noise_array,MRSCont.files{kk},kk);
end
%% perform SENSE unfolding
for kk = 1:MRSCont.nDatasets
    MRSCont.raw{kk} = osp_SenseUnfolding(MRSCont.raw_uncomb{kk},MRSCont.SENSE{kk});
    if MRSCont.flags.hasRef
        MRSCont.raw_ref{kk} = osp_SenseUnfolding(MRSCont.raw_ref_uncomb{kk},MRSCont.SENSE{kk});
    end
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        MRSCont.raw_w{kk} = osp_SenseUnfolding(MRSCont.raw_w_uncomb{kk},MRSCont.SENSE{kk});
    end
end

% Set flags
MRSCont.flags.coilsCombined     = 1;
% Delete un-combined data to free up memory
raw_fields = {'raw_uncomb','raw_ref_uncomb','raw_w_uncomb'};
for kk = 1:length(raw_fields)
    if isfield(MRSCont, raw_fields{kk})
        MRSCont = rmfield(MRSCont, raw_fields{kk});
    end
end
end
