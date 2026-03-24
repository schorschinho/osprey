function  [niiSOS, niiSENSE, noise_array] = osp_cpx2nii(cpxFile, sinFile, niiSOSFile, niiSENSEFile, launchViewer)
%%[niiSOS, niiSENSE] = osp_cpx2nii(cpxFile, sinFile, niiSOSFile, niiSENSEFile, launchViewer)
%   Reads complex Philips data (*.cpx) and geometry (.sin) files and
%   converts the image to NIfTI (*.nii) format.
%
%   Input:
%       cpxFile = Complex Philips raw image data file
%       sinFile = Accompanying Philips sin file
%       niiSOSFile = Output file name under which the SOS file will be saved
%       niiSENSEFile = Output file name under which the 4D coil image will
%       be saved
%       launchViewer = Start external nii viewer to view output?
%
%   Output:
%       niiSOS = 3D NIfTI containing the sum-of-squares image
%       niiSENSE = 4D NIfTI containing the complex images from each coil
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2020-12-21)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This function requires the DICOM to NIfTI conversion toolbox
%       https://github.com/xiangruili/dicm2nii
%       (Xiangrui Li, Ohio State University)
%       to be added to the MATLAB path.
%
%   History:
%       2020-12-21: First version of the code.
%

if nargin < 5
    launchViewer = 0;
end

% Extract geometry information
cpxSinInfo    = get_sin_info(sinFile);

% Load complex image data
[cpx_file, Ref, img_array, img_body, noise_array, noise_body] = loadRefScan(cpxFile);

% Extract ref scan orientation
refsc_orientation = cpxSinInfo.slice_orientation;
refsc_voxel_size  = cpxSinInfo.voxel_size;
refsc_offset      = cpxSinInfo.loc_ap_rl_fh_offcentres;

% Calculate sum of squares
SOS = squeeze(sqrt(sum(img_array.*conj(img_array),1)));
% Generate sensitivity maps
% Now create the individual complex coil sensitivity maps by normalizing
% every individual coil image to the sum-of-squares image
% Note: Outside of the brain, there is going to be noise in the SOS image,
% so we'll zero that out in order to get better behaved sensitivity maps.
sensitivities = zeros(size(img_array));
mask = (SOS/max(SOS(:)))>0.04;
for c=1:size(sensitivities,1)
    sensitivities(c,:,:,:) = mask.*squeeze(img_array(c,:,:,:))./SOS;
end

% We know the FoV
FoV(1) = size(SOS,1)*refsc_voxel_size(1);
FoV(2) = size(SOS,2)*refsc_voxel_size(2);
FoV(3) = size(SOS,3)*refsc_voxel_size(3);

switch refsc_orientation
    case 2 % coronal slices
        % This arrives in the order
        % [ncoils ap fh rl]
        SOS = permute(SOS, [2 3 1]);
        sensitivities = permute(sensitivities, [3 4 2 1]);
        
        % Voxel sizes are not in AP-RL-FH order, need to be permuted
        % But stack offset comes in that order, does not need permutation
        FoV = [FoV(2) FoV(3) FoV(1)];
        refsc_voxel_size = [refsc_voxel_size(2) refsc_voxel_size(3) refsc_voxel_size(1)];
        
    case 1 % sagittal slices
        
        % This arrives in the order
        % [ncoils fh ap rl]
        SOS = permute(SOS,[3 2 1]);
        sensitivities = permute(sensitivities, [4 3 2 1]);
        
        % Voxel sizes are not in AP-RL-FH order, need to be permuted
        % But stack offset comes in that order, does not need permutation
        FoV = [FoV(3) FoV(2) FoV(1)];
        refsc_voxel_size = [refsc_voxel_size(3) refsc_voxel_size(2) refsc_voxel_size(1)];
        
    case 0 % axial slices
        
        % Essentially, do nothing!
        
        % This arrives in the order
        % [ncoils rl ap fh]
        SOS = permute(SOS,[1 2 3]);
        sensitivities = permute(sensitivities, [2 3 4 1]);
        
        % Voxel sizes are not in AP-RL-FH order, need to be permuted
        % But stack offset comes in that order, does not need permutation
        FoV = [FoV(1) FoV(2) FoV(3)];
        refsc_voxel_size = [refsc_voxel_size(1) refsc_voxel_size(2) refsc_voxel_size(3)];

end

% But we need to find the origin of the stack from the center coordinates
% We know the center of the stack
stackCenter = refsc_offset;

% Calculate corners of the FOV
FoV_corners = ...
    [FoV(1)/2 -FoV(2)/2  FoV(3)/2;
    -FoV(1)/2 -FoV(2)/2  FoV(3)/2;
    -FoV(1)/2  FoV(2)/2  FoV(3)/2;
     FoV(1)/2  FoV(2)/2  FoV(3)/2;
    -FoV(1)/2  FoV(2)/2 -FoV(3)/2;
     FoV(1)/2  FoV(2)/2 -FoV(3)/2;
     FoV(1)/2 -FoV(2)/2 -FoV(3)/2;
    -FoV(1)/2 -FoV(2)/2 -FoV(3)/2];

FoV_corners = FoV_corners + stackCenter;

% Now we build the nifti array
% Flip LR and AP directions
SOS = flip(SOS, 1);
SOS = flip(SOS, 2);
sensitivities = flip(sensitivities, 1);
sensitivities = flip(sensitivities, 2);

% Initialize the nii header container including the header
% Fill in the necessary geometry fields in the header
% Pixel dimensions (LR - AP - FH)
% First save the sum-of-squares only
niiSOS = nii_tool('init', SOS);
niiSOS.hdr.pixdim = [1 refsc_voxel_size(1) refsc_voxel_size(2) refsc_voxel_size(3) 1 1 1 1];
%Then save the sensitivities as 4D array (with the number of coils in
% the 4-th dimension, analogous to fMRI time series)
niiSENSE = nii_tool('init', sensitivities);
niiSENSE.hdr.pixdim = [1 refsc_voxel_size(1) refsc_voxel_size(2) refsc_voxel_size(3) 0 1 1 1];

%nii.hdr.xyzt_units = 2;
%nii.hdr.qform_code = 1;
%nii.hdr.quatern_b = 0;
%nii.hdr.quatern_c = 0;
%nii.hdr.quatern_d = -1;
%nii.hdr.qoffset_x = refsc_offset(1);
%nii.hdr.qoffset_y = refsc_offset(2);
%nii.hdr.qoffset_z = refsc_offset(3);
% below is work in progress
r = 6;
niiSOS.hdr.sform_code = 2;
niiSOS.hdr.srow_x = [refsc_voxel_size(1) 0 0 -FoV_corners(r,1)+refsc_voxel_size(1)/2];
niiSOS.hdr.srow_y = [0 refsc_voxel_size(2) 0 -FoV_corners(r,2)+refsc_voxel_size(2)/2];
niiSOS.hdr.srow_z = [0 0 refsc_voxel_size(3) FoV_corners(r,3)+refsc_voxel_size(3)/2];
niiSOS.hdr.intent_name = '';
niiSOS.hdr.magic = 'n+1';
niiSOS.hdr.version = 1;
niiSOS.hdr.swap_endian = 0;
niiSENSE.hdr.sform_code = 2;
niiSENSE.hdr.srow_x = [refsc_voxel_size(1) 0 0 -FoV_corners(r,1)+refsc_voxel_size(1)/2];
niiSENSE.hdr.srow_y = [0 refsc_voxel_size(2) 0 -FoV_corners(r,2)+refsc_voxel_size(2)/2];
niiSENSE.hdr.srow_z = [0 0 refsc_voxel_size(3) FoV_corners(r,3)+refsc_voxel_size(3)/2];
niiSENSE.hdr.intent_name = '';
niiSENSE.hdr.magic = 'n+1';
niiSENSE.hdr.version = 1;
niiSENSE.hdr.swap_endian = 0;



nii_tool('save', niiSOS, niiSOSFile);
nii_tool('save', niiSENSE, niiSENSEFile);

% Optionally: launch the external viewing tool
if launchViewer
    nii_viewer(niiSOSFile);
    nii_viewer(niiSENSEFile);
end

end
