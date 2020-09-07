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
[~] = osp_Toolbox_Check('OspreyCoreg',MRSCont.flags.isGUI);
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
                        msg = 'Philips DATA files do not contain voxel geometry information.';
                        fprintf(fileID,msg);
                        error(msg);                        
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


% If a pair of secondary images exists, start the process of
% co-registering the first T1 to the second T1 here.
% We choose to keep the pair of secondary images in native space,
% because these are likely of lower resolution (for example PET or DWI
% data), and we don't want to reslice them.
for kk = 1:MRSCont.nDatasets
    if (MRSCont.flags.hasSecondT1 && MRSCont.flags.hasPET)
            
        % Unfortunately, we need to temporarily duplicate the initial T1 we
        % want to work on, since the SPM functions change the header!
        originalT1File = MRSCont.coreg.vol_image{kk}.fname;
        [originalT1Path, originalT1Filename, originalT1Ending] = fileparts(originalT1File);
        originalT1Copy = fullfile(originalT1Path, [originalT1Filename '_tempCopy' originalT1Ending]);
        copyfile(originalT1File, originalT1Copy);
        % Also copy the voxel mask
        originalVoxelMaskFile = MRSCont.coreg.vol_mask{kk}.fname;
        [originalVoxelMaskPath, originalVoxelMaskFilename, originalVoxelMaskEnding] = fileparts(originalVoxelMaskFile);
        originalVoxelMaskCopy = fullfile(originalVoxelMaskPath, [originalVoxelMaskFilename '_tempCopy' originalVoxelMaskEnding]);
        copyfile(originalVoxelMaskFile, originalVoxelMaskCopy);
        
        % Set up the SPM12 coregister and re-slice batch
        % Choose the second T1 provided in the job file as target (i.e.
        % stationary image, since this is already co-registered with the
        % secondary modality image (e.g. PET) in files_pet.
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {MRSCont.files_nii2{kk}};
        % Choose the original T1 copy we just created as the 'moving'
        % image.
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {originalT1Copy};
        % We will also apply the rigid transformation to the copy of the 
        % voxel mask we just created.
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {originalVoxelMaskCopy};
        % Set up default options for the estimation and reslicing.
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        % Run the job!
        spm_jobman('run',matlabbatch);
        
        % Get filenames of the voxel mask resliced to the secondary T1 -
        % this is easy, SPM just prefixed it with an 'r'.
        % Rename it so that it does not contain the '_tempCopy' bit any
        % more.
        reslicedVoxelMaskCopy = fullfile(originalVoxelMaskPath, ['r' originalVoxelMaskFilename '_tempCopy' originalVoxelMaskEnding]);
        reslicedVoxelMaskFile = fullfile(originalVoxelMaskPath, ['r' originalVoxelMaskFilename originalVoxelMaskEnding]);
        movefile(reslicedVoxelMaskCopy, reslicedVoxelMaskFile)
        
        % Create SPM volume and read in the NIfTI file with the secondary T1.
        vol_image_2nd   = spm_vol(MRSCont.files_nii2{kk});
        vol_image_pet   = spm_vol(MRSCont.files_pet{kk});
        vol_mask_2nd    = spm_vol(reslicedVoxelMaskFile);
        [T1_2nd,XYZ]    = spm_read_vols(vol_image_2nd);
        [PET,XYZ]       = spm_read_vols(vol_image_pet);
        T1_2nd_max      = max(T1_2nd(:));
        PETmax          = max(PET(:));
        MRSCont.coreg.vol_image_2nd{kk} = vol_image_2nd;
        MRSCont.coreg.T1_2nd_max{kk}    = T1_2nd_max;
        MRSCont.coreg.vol_image_pet{kk} = vol_image_pet;
        MRSCont.coreg.PETmax{kk}        = PETmax;
        MRSCont.coreg.vol_mask_2nd{kk}  = vol_mask_2nd;
        
        %%% 3. DETERMINE PET IMAGE INTENSITY %%%
        % Get the histogram of PET image intensities inside the voxel
        % Get the input file name
        [path,filename,~]   = fileparts(MRSCont.files{kk});
        [maskDir, maskName, maskExt] = fileparts(vol_mask_2nd.fname);
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_split          = regexp(path,filesep,'split');
        if length(path_split) > 2
            saveName = [path_split{end-1} '_' path_split{end} '_' filename];
        end
        vol_PETMask.fname    = fullfile(saveDestination, [saveName '_PET' maskExt]);
        vol_PETMask.descrip  = 'PETmasked_MRS_Voxel_Mask';
        vol_PETMask.dim      = vol_mask_2nd.dim;
        vol_PETMask.dt       = vol_mask_2nd.dt;
        vol_PETMask.mat      = vol_mask_2nd.mat;
        PET_voxmask_vol      = vol_image_pet.private.dat(:,:,:) .* vol_mask_2nd.private.dat(:,:,:);
        vol_PETMask          = spm_write_vol(vol_PETMask, PET_voxmask_vol);
        
        % Get everything greater than zero
        PETIntensityInVoxel = nonzeros(vol_PETMask.private.dat(:,:,:));
        rawPETIntensitySum = sum(PETIntensityInVoxel);
        
        idx = PETIntensityInVoxel > 0.2;
        PETIntensityInVoxel = PETIntensityInVoxel(idx);
        nBins = round(sqrt(length(PETIntensityInVoxel)));           % Number of bins
        [yIntensity, xIntensity] = hist(PETIntensityInVoxel,nBins);
        
        % Initial guesses
        [maxIntensityInit, maxIndexInit] = max(yIntensity);
        % Get starting values for linear baseline and offset
        grad_points = (yIntensity(end) - yIntensity(1)) ./ abs(xIntensity(end) - xIntensity(1));
        LinearInit = grad_points ./ abs(xIntensity(1) - xIntensity(2));
        constInit = (yIntensity(end)+yIntensity(1))/2;
        
        % Set initials and run model
        GaussModelInit = [maxIntensityInit -0.2 xIntensity(maxIndexInit) -LinearInit constInit];
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6,'FunValCheck','off');
        [GaussModelParams, residPlot] = nlinfit(xIntensity, yIntensity, @GaussModel, GaussModelInit, nlinopts); % re-run for residuals for output figure

        
        % Save all kinds of estimates
        MRSCont.coreg.pet.rawPETIntensitySum(kk) = rawPETIntensitySum;
        MRSCont.coreg.pet.histogram.xIntensity{kk} = xIntensity;
        MRSCont.coreg.pet.histogram.yIntensity{kk} = yIntensity;
        MRSCont.coreg.pet.histogram.fitParams{kk}  = GaussModelParams;
        MRSCont.coreg.pet.histogram.mostFrequentIntensity(kk)  = GaussModelParams(3);
        MRSCont.coreg.pet.histogram.distFWHM(kk)  = GaussModelParams(2);
        
        % Delete the temporary copy of the voxel masks
        delete(sprintf(originalVoxelMaskCopy));
        
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

% Delete the temporary copies and reslices of the T1
delete(sprintf(originalT1Copy));
reslicedT1File = fullfile(originalT1Path, ['r' originalT1Filename '_tempCopy' originalT1Ending]);
delete(sprintf(reslicedT1File));

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