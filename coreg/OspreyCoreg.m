function [MRSCont] = OspreyCoreg(MRSCont)
%% [MRSCont] = OspreyCoreg(MRSCont)
%   This function parses header information about the dimensions, location, 
%   and rotations of the MRS voxel. It then proceeds to call SPM12 to 
%   import structural images that have been provided in the job file.
%   These images are supplied in NIfTI format (*.nii) for Philips and
%   Siemens, or as DICOM folders for GE data.
%
%   OspreyCoreg then uses the voxel geometry information to create a
%   voxel mask (in NIfTI format).
%
%   USAGE:
%       MRSCont = OspreyCoreg(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-28)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-28: First version of the code.

% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to process data, but raw data has not been loaded yet. Run OspreyLoad first.')
end

% Version check
MRSCont.ver.CheckCoreg            = '100 Coreg';

warning('off','all');

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'VoxelMasks');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

%% Loop over all datasets
refProcessTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressbar = waitbar(0,'Start','Name','Osprey Coregister');
    waitbar(0,progressbar,sprintf('Coregistered voxel from dataset %d out of %d total datasets...\n', 0, MRSCont.nDatasets))
end
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Coregistering voxel from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    if ((MRSCont.flags.didCoreg == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'coreg') && (kk > length(MRSCont.coreg.vol_image))) || ~isfield(MRSCont.ver, 'Coreg') || ~strcmp(MRSCont.ver.Coreg,MRSCont.ver.CheckCoreg))

        % Get the input file name
        [path,filename,~]   = fileparts(MRSCont.files{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_split          = regexp(path,filesep,'split');
        if length(path_split) > 2
            saveName = [path_split{end-1} '_' path_split{end} '_' filename];
        end
        % Generate file name for the voxel mask NIfTI file to be saved under
        maskFile            = fullfile(saveDestination, [saveName '_VoxelMask.nii']);

        % Call voxel mask generator depending on file type
        switch MRSCont.vendor
            case 'Siemens'
                % Load the *.nii file provided in the job file
                vol_image = spm_vol(MRSCont.files_nii{kk});
                switch MRSCont.datatype
                    case 'TWIX'
                        [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{kk}, vol_image, maskFile);
                    case 'RDA'
                        [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{kk}, vol_image, maskFile);
                    case 'DICOM'
                        [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{kk}, vol_image, maskFile);
                    otherwise
                        error('Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).');
                end
            case 'Philips'
                % Load the *.nii file provided in the job file
                vol_image = spm_vol(MRSCont.files_nii{kk});
                switch MRSCont.datatype
                    case 'SDAT'
                        [vol_mask, T1_max, voxel_ctr] = coreg_sdat(MRSCont.raw{kk}, vol_image, maskFile);
                    case 'DATA'
                        error('Philips DATA files do not contain voxel geometry information.');
                    case 'RAW'
                        error('Philips RAW files do not contain voxel geometry information.');
                    otherwise
                        error('Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).');
                end
            case 'GE'
                switch MRSCont.datatype
                    case 'P'
                        % Load the DICOM folder provided in the job file
                        dcm_folder = MRSCont.files_nii{kk};
                        [vol_mask, T1_max, vol_image, voxel_ctr] = coreg_p(MRSCont.raw{kk}, dcm_folder, maskFile);
                    otherwise
                        error('Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).');
                end
            otherwise
                error('Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).');
        end

        % Save back the image and voxel mask volumes to MRSCont
        MRSCont.coreg.vol_image{kk} = vol_image;
        MRSCont.coreg.vol_mask{kk}  = vol_mask;
        MRSCont.coreg.T1_max{kk}    = T1_max;
        MRSCont.coreg.voxel_ctr{kk} = voxel_ctr;
    end
    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Coregistered voxel from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets))
    end
end
fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    close(progressbar)
end
toc(refProcessTime);

%% Clean up and save
% Set exit flags and version
MRSCont.flags.didCoreg           = 1;
MRSCont.ver.Coreg            = '100 Coreg';

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if ~MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont');
end

end