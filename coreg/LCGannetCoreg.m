function [MRSCont] = LCGannetCoreg(MRSCont)
%% [MRSCont] = LCGannetCoreg(MRSCont)
%   This function calls SPM12 to handle an image provided in NIfTI format
%   (*.nii), and then uses the voxel geometry information provided in the
%   MRS data headers to create a voxel mask (also in NIfTI format).
%
%   USAGE:
%       MRSCont = LCGannetCoreg(MRSCont);
%
%   INPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
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

% Close any remaining open figures
close all;
warning('off','all');

%% Loop over all datasets
refProcessTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Coregistering voxel from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    % Load the *.nii file provided in LCG_filesToLoad
    vol_image = spm_vol(MRSCont.files_nii{kk});
    
    % Call voxel mask generator depending on file type
    switch MRSCont.vendor
        case 'Siemens'
            switch MRSCont.datatype
                case 'TWIX'
                    [vol_mask, T1_max] = coreg_siemens(MRSCont.raw{kk}, vol_image);
                case 'RDA'
                    [vol_mask, T1_max] = coreg_siemens(MRSCont.raw{kk}, vol_image);
                case 'DICOM'
                    [vol_mask, T1_max] = coreg_siemens(MRSCont.raw{kk}, vol_image);
                otherwise
                    error('Data type not supported. Please contact the LCGannet team (gabamrs@gmail.com).');
            end
        case 'Philips'
            switch MRSCont.datatype
                case 'SDAT'
                    [vol_mask, T1_max] = coreg_sdat(MRSCont.raw{kk}, vol_image);
                case 'DATA'
                    error('Philips DATA files do not contain voxel geometry information. Please enter the corresponding SPAR files in LCG_filesToLoad.');
                case 'RAW'
                    error('Philips RAW files do not contain voxel geometry information. Please enter the corresponding SPAR files in LCG_filesToLoad.');
                otherwise
                    error('Data type not supported. Please contact the LCGannet team (gabamrs@gmail.com).');
            end
        case 'GE'
            switch MRSCont.datatype
                case 'P'
                    [vol_mask, T1_max] = coreg_p(MRSCont.raw{kk}, vol_image);
                otherwise
                    error('Data type not supported. Please contact the LCGannet team (gabamrs@gmail.com).');
            end
        otherwise
            error('Vendor not supported. Please contact the LCGannet team (gabamrs@gmail.com).');
    end
    
    % Save back the image and voxel mask volumes to MRSCont
    MRSCont.coreg.vol_image{kk} = vol_image;
    MRSCont.coreg.vol_mask{kk}  = vol_mask;
    MRSCont.coreg.T1_max{kk}    = T1_max;
end
fprintf('... done.\n');
toc(refProcessTime);

%% Clean up and save
% Set exit flags
MRSCont.flags.didCoreg           = 1;

end