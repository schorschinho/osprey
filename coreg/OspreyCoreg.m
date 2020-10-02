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

outputFolder = MRSCont.outputFolder;
fileID = fopen(fullfile(outputFolder, 'LogFile.txt'),'a+');
% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    msg = 'Trying to process data, but raw data has not been loaded yet. Run OspreyLoad first.';
    fprintf(fileID,msg);
    error(msg);
end

% Version, toolbox check and updating log file
MRSCont.ver.CheckCoreg       = '1.0.0 Coreg';
fprintf(fileID,['Timestamp %s ' MRSCont.ver.Osp '  ' MRSCont.ver.CheckCoreg '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
[~] = osp_Toolbox_Check ('OspreyCoreg',MRSCont.flags.isGUI);

warning('off','all');

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'VoxelMasks');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

%% Loop over all datasets
refCoregTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Coregistering voxel from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    fprintf(fileID,[reverseStr, msg]);
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('Coregistering voxel from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
        drawnow
    end
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
                        msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                        fprintf(fileID,msg);
                        error(msg);
                end
            case 'Philips'
                % Load the *.nii file provided in the job file
                vol_image = spm_vol(MRSCont.files_nii{kk});
                switch MRSCont.datatype
                    case 'SDAT'
                        [vol_mask, T1_max, voxel_ctr] = coreg_sdat(MRSCont.raw{kk}, vol_image, maskFile);
                    case 'DATA'
                        if isfield(MRSCont.raw{kk}, 'geometry')
                            [vol_mask, T1_max, voxel_ctr] = coreg_sdat(MRSCont.raw{kk}, vol_image, maskFile);
                        else
                        msg = 'Philips DATA files do not contain voxel geometry information.';
                        fprintf(fileID,msg);
                        error(msg);  
                      end
                    case 'RAW'
                        msg = 'Philips RAW files do not contain voxel geometry information.';
                        fprintf(fileID,msg);
                        error(msg);                        
                    otherwise
                        msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                        fprintf(fileID,msg);
                        error(msg);                        
                end
            case 'GE'
                [~,~,file_exten]=fileparts(MRSCont.files_nii{kk});
                if contains(file_exten,'.nii')
                    % Load the *.nii file provided in the job file
                    vol_image = spm_vol(MRSCont.files_nii{kk});
                    
                    [vol_mask, T1_max, voxel_ctr] = coreg_ge_nifti(MRSCont.raw{kk}, vol_image, maskFile);
                else
                    switch MRSCont.datatype
                        case 'P'
                            % Load the DICOM folder provided in the job file
                            dcm_folder = MRSCont.files_nii{kk};
                            [vol_mask, T1_max, vol_image, voxel_ctr] = coreg_p(MRSCont.raw{kk}, dcm_folder, maskFile);
                        otherwise
                            msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                            fprintf(fileID,msg);
                            error(msg);  
                    end
                end
            otherwise
                msg = 'Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(fileID,msg);
                error(msg);                
        end

        % Save back the image and voxel mask volumes to MRSCont
        MRSCont.coreg.vol_image{kk} = vol_image;
        MRSCont.coreg.vol_mask{kk}  = vol_mask;
        MRSCont.coreg.T1_max{kk}    = T1_max;
        MRSCont.coreg.voxel_ctr{kk} = voxel_ctr;
    end
end
fprintf('... done.\n');
time = toc(refCoregTime);
if MRSCont.flags.isGUI        
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf(fileID,'... done.\n Elapsed time %f seconds\n',time);
MRSCont.runtime.Coreg = time;
fclose(fileID); %close log file
%% Clean up and save
% Set exit flags and version
MRSCont.flags.didCoreg           = 1;
MRSCont.ver.Coreg            = '1.0.0 Coreg';

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    for kk = 1 : MRSCont.nDatasets
            osp_plotModule(MRSCont, 'OspreyCoreg', kk);
    end
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end