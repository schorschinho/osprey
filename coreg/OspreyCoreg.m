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
diary(fullfile(outputFolder, 'LogFile.txt'));

warning('off','all');
% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyCoreg');

%Do some check on the naming convention first to avoid overwriting the
%output due to non BIDS conform data names
SameName = 0;
if MRSCont.nDatasets(1) > 1
    specFile = MRSCont.files{1,1};
    specFile2 = MRSCont.files{1,2};
    [~, SpecName, ~]  = fileparts(specFile);
    [~, SpecName2, ~]  = fileparts(specFile2);
    if ~isempty(SpecName) && ~isempty(SpecName2)
        SameName = strcmp(SpecName,SpecName2);
    else
        [DirName, ~, ~]  = fileparts(specFile);
        [DirName2, ~, ~]  = fileparts(specFile2);
        SepFiles =  split(DirName, filesep);
        SepFiles(strcmp(SepFiles,''))=[];
        DirName = SepFiles{end};
        SepFiles =  split(DirName2, filesep);
        SepFiles(strcmp(SepFiles,''))=[];
        DirName2 = SepFiles{end};
        SameName = strcmp(DirName,DirName2);
    end
end
MRSCont.coreg.SameName = SameName;

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'VoxelMasks'); %CWDJ - Address in future update
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

%% Loop over all datasets
refCoregTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
for kk = 1:MRSCont.nDatasets(1)
     [~] = printLog('OspreyCoreg',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);  
    if ~(MRSCont.flags.didCoreg == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'coreg') && (kk > length(MRSCont.coreg.vol_image))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
        
        if SameName
            [PreFix] = osp_generate_SubjectAndSessionPrefix(MRSCont.files{1,kk},kk);
            PreFix = [PreFix '_'];
        else
            PreFix = '';
        end

        % Get the input file name
        if ~exist('DirName','var')
            [~,filename,fileext]   = fileparts(MRSCont.files{1,kk});
            if strcmp(fileext,'.gz')
                [~,filename,~]   = fileparts(filename);
            end
        else
            [dirname,~,~]   = fileparts(MRSCont.files{1,kk});
            SepFiles =  split(dirname, filesep);
            SepFiles(strcmp(SepFiles,''))=[];
            filename = SepFiles{end};
        end
        % Get the nii file name
        [~,~,T1ext]   = fileparts(MRSCont.files_nii{kk});
        
        %<source_entities>[_space-<space>][_res-<label>][_den-<label>][_label-<label>][_desc-<label>]_mask.nii.gz
        saveName = [PreFix filename '_space-scanner']; %CWDJ Check space.
        
        % Generate file name for the voxel mask NIfTI file to be saved under
        maskFile            = fullfile(saveDestination, [saveName '_mask.nii']);

        %Uncompress .nii.gz if needed
        if strcmp(T1ext,'.gz')
            gunzip(MRSCont.files_nii{kk});
            MRSCont.files_nii{kk} = strrep(MRSCont.files_nii{kk},'.gz','');
        end

        if strcmp(MRSCont.datatype,'NIfTI-MRS')
            vol_image = spm_vol(MRSCont.files_nii{kk});
            [vol_mask, T1_max, voxel_ctr] = coreg_nifti(MRSCont.raw{1,kk}, vol_image, maskFile);
        else
            % Call voxel mask generator depending on file type
            switch MRSCont.vendor
                case 'Siemens'
                    % Load the *.nii file provided in the job file
                    
                    vol_image = spm_vol(MRSCont.files_nii{kk});
                    switch MRSCont.datatype
                        case 'TWIX'
                            [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{1,kk}, vol_image, maskFile);
                        case 'RDA'
                            [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{1,kk}, vol_image, maskFile);
                        case 'DICOM'
                            [vol_mask, T1_max, voxel_ctr] = coreg_siemens(MRSCont.raw{1,kk}, vol_image, maskFile);
                        otherwise
                            msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                            fprintf(msg);
                            error(msg);
                    end
                case 'Philips'
                    % Load the *.nii file provided in the job file
                    vol_image = spm_vol(MRSCont.files_nii{kk});
                    switch MRSCont.datatype
                        case 'SDAT'
                             if ~MRSCont.flags.isMRSI % SVS coregistration
                                [vol_mask, T1_max, voxel_ctr,~] = coreg_sdat(MRSCont.raw{1,kk}, vol_image, maskFile);
                             else
                                  [vol_mask, T1_max, voxel_ctr,~] = coreg_sdat(MRSCont.raw{1,kk}, vol_image, maskFile,2);
    %                                MRSCont.coreg.vol_mask_mrsi{kk} = vol_mask_mrsi;
                             end
                        case 'DATA'
                            if isfield(MRSCont.raw{kk}, 'geometry')
                                if ~MRSCont.flags.isPRIAM % SVS coregistration
                                    [vol_mask, T1_max, voxel_ctr,~] = coreg_sdat(MRSCont.raw{1,kk}, vol_image, maskFile);
                                else 
                                    [vol_mask, T1_max, voxel_ctr,~] = coreg_sdat(MRSCont.raw{1,kk}, vol_image, maskFile, MRSCont.SENSE{kk});
                                end
                            else
                            msg = 'Philips DATA files do not contain voxel geometry information.';
                            fprintf(msg);
                            error(msg);  
                            end
                        case 'RAW'
                            msg = 'Philips RAW files do not contain voxel geometry information.';
                            fprintf(msg);
                            error(msg);                        
                        otherwise
                            msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                            fprintf(msg);
                            error(msg);                        
                    end
                case 'GE'
                    [~,~,file_exten]=fileparts(MRSCont.files_nii{kk});
                    if contains(file_exten,'.nii')
                        % Load the *.nii file provided in the job file
                        vol_image = spm_vol(MRSCont.files_nii{kk});
                        
                        [vol_mask, T1_max, voxel_ctr] = coreg_ge_nifti(MRSCont.raw{1,kk}, vol_image, maskFile);
                    else
                        switch MRSCont.datatype
                            case 'P'
                                % Load the DICOM folder provided in the job file                           
                                [vol_mask, T1_max, vol_image, voxel_ctr] = coreg_p(MRSCont.raw{1,kk}, MRSCont.files_nii{kk}, maskFile);
                            otherwise
                                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                                fprintf(msg);
                                error(msg);  
                        end
                    end
                otherwise
                    msg = 'Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                    fprintf(msg);
                    error(msg);                
            end
        end

        if MRSCont.opts.img.deface
            [anaon_vol_image] = spm_deface({vol_image.fname});
            saveDestinationAnon = fullfile(MRSCont.outputFolder, 'DefacedNII'); %CWDJ - Address in future update
            if ~exist(saveDestinationAnon,'dir')
                mkdir(saveDestinationAnon);
            end
            [~,AnonName,AnonExt] = fileparts(anaon_vol_image{1});
            [~] = movefile(anaon_vol_image{1}, fullfile(saveDestinationAnon,[AnonName AnonExt]));
            vol_image = spm_vol(fullfile(saveDestinationAnon,[AnonName AnonExt]));
        end

        % Save back the image and voxel mask volumes to MRSCont
        MRSCont.coreg.vol_image{kk} = vol_image;
        MRSCont.coreg.vol_mask{kk}  = vol_mask;
        MRSCont.coreg.T1_max{kk}    = T1_max;
        MRSCont.coreg.voxel_ctr{kk} = voxel_ctr;
        
        if MRSCont.flags.addImages
            [MRSCont.coreg.three_plane_img{kk}] = osp_extract_three_plane_image(vol_image, vol_mask,voxel_ctr,T1_max);
        end
        
        %Delete .nii file if a .nii.gz
         if strcmp(T1ext,'.gz')
            delete(MRSCont.files_nii{kk});
            MRSCont.files_nii{kk} = strrep(MRSCont.files_nii{kk},'.nii','.nii.gz');
         end
         if ~MRSCont.flags.isPRIAM
             gzip(vol_mask.fname);
             delete(vol_mask.fname);
         else
             gzip(vol_mask{1}.fname);
             delete(vol_mask{1}.fname);
             gzip(vol_mask{2}.fname);
             delete(vol_mask{2}.fname);
         end
            
    end
end
time = toc(refCoregTime);
[~] = printLog('done',time,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
MRSCont.runtime.Coreg = time;
%% Clean up and save
% Set exit flags and version
MRSCont.flags.didCoreg           = 1;
diary off

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyCoreg')
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end