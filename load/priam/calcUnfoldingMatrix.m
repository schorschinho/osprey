function MRSCont = calcUnfoldingMatrix(MRSCont, niiSOSFile, niiSENSEFile, kk, launchViewer)
%% function MRSCont = calcUnfoldingMatrix(MRSCont, niiSOSFile, niiSENSEFile, kk, launchViewer)
%   Calculates the unfolding matrix for a given coil reference image and
%   voxel locations specified in MRS_struct.
%
%   Input:
%       MRSCont     = Osprey MRS data container.
%
%   Output:
%       MRSCont     = Osprey MRS data container filled with information 
%                   about the SENSE reconstruction (in MRSCont.SENSE)
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-18)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-18: First version of the code.
%

if nargin < 4
    launchViewer = 0;
end


% Deactivate MATLAB warnings and load geometry parameters
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');


%%% 1. GENERATE THE COORDINATES OF THE VOXEL CORNERS
% Get dimensions, orientation and global offset from MRS scan
% First, determine the filename and make sure the *.sin file exists,
% otherwise throw an error
sin_file = strrep(MRSCont.files{kk},'.data','.sin');
if ~isfile(sin_file)
    msg = 'No *.sin file for the MRS scan found. It needs to have the same name as the .data file.';
    error(msg);
end
% Now actually extract the information we need
MRSCont.SENSE{kk}.sin_info  = get_sin_info(sin_file);
% Get geometry info from SIN and SDAT files
lr_size = MRSCont.SENSE{kk}.sin_info.voxel_size(2); % saved in sin as row - column - slice
ap_size = MRSCont.SENSE{kk}.sin_info.voxel_size(1); 
cc_size = MRSCont.SENSE{kk}.sin_info.voxel_size(3);
lr_off  = MRSCont.SENSE{kk}.sin_info.loc_ap_rl_fh_offcentres(2); % saved in sin as ap - rl - fh
ap_off  = MRSCont.SENSE{kk}.sin_info.loc_ap_rl_fh_offcentres(1); 
cc_off  = MRSCont.SENSE{kk}.sin_info.loc_ap_rl_fh_offcentres(3);
lr_ang  = MRSCont.raw_uncomb{kk}.geometry.rot.lr; % not saved in sin at all! So need to pull in from the loaded SDAT information header 
ap_ang  = MRSCont.raw_uncomb{kk}.geometry.rot.ap; 
cc_ang  = MRSCont.raw_uncomb{kk}.geometry.rot.cc;

% We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;
ap_ang = -ap_ang;
lr_ang = -lr_ang;

% Define the voxel - use x y z
% Currently have spar convention that have in AUD voxel - will need to
% check for everything in future...
% x - left = positive
% y - posterior = postive
% z - superior = positive
vox_ctr = ...
    [lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2  cc_size/2;
     lr_size/2  ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2 -ap_size/2 -cc_size/2;
    -lr_size/2 -ap_size/2 -cc_size/2];

% Make rotations on voxel
rad = pi/180;
initrot = zeros(3,3);

xrot      = initrot;
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang *rad);
xrot(2,3) = -sin(lr_ang*rad);
xrot(3,2) = sin(lr_ang*rad);
xrot(3,3) = cos(lr_ang*rad);

yrot      = initrot;
yrot(1,1) = cos(ap_ang*rad);
yrot(1,3) = sin(ap_ang*rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang*rad);
yrot(3,3) = cos(ap_ang*rad);

zrot      = initrot;
zrot(1,1) = cos(cc_ang*rad);
zrot(1,2) = -sin(cc_ang*rad);
zrot(2,1) = sin(cc_ang*rad);
zrot(2,2) = cos(cc_ang*rad);
zrot(3,3) = 1;

% Rotate voxel
vox_rot = xrot * yrot * zrot * vox_ctr.';

% Calculate corner coordinates in real space relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;

% Calculate unit vectors pointing along the voxel edges in real space
lr_vec = vox_corner(:,1) - vox_corner(:,2);
ap_vec = vox_corner(:,4) - vox_corner(:,1);
cc_vec = vox_corner(:,1) - vox_corner(:,7);
lr_vec = lr_vec ./ (norm(lr_vec));
ap_vec = ap_vec ./ (norm(ap_vec));
cc_vec = cc_vec ./ (norm(cc_vec));

% Depending on the PRIAM offset direction parameter, calculate the voxel
% corner and voxel center coordinates
secondVoxelOffset    = MRSCont.SENSE{kk}.priam_offset;
secondVoxelDirection = MRSCont.SENSE{kk}.priam_direction;
switch secondVoxelDirection
    case 'R'
        shift = secondVoxelOffset*lr_vec;
    case 'L'
        shift = -secondVoxelOffset*lr_vec;
    case 'A'
        shift = secondVoxelOffset*ap_vec;
    case 'P'
        shift = -secondVoxelOffset*ap_vec;
    case 'F'
        shift = -secondVoxelOffset*cc_vec;
    case 'H'
        shift = -secondVoxelOffset*cc_vec;
end
vox_corner(:,:,2) = vox_corner + shift;
vox_center = [lr_off ap_off cc_off]';
vox_center(:,2) = vox_center + shift;


%%% 2. PREPARE THE COIL SENSITIVITY IMAGE
% Create SPM volume and read in the NIfTI file with the coil sensitivity sum-of-squares image.
vol_sos      = spm_vol(niiSOSFile);
[SOS,XYZ]    = spm_read_vols(vol_sos);
SOS_max      = max(SOS(:));
%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(vol_sos,'fv');
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);


%%% 3. CREATE AND SAVE THE VOXEL MASKS
% Loop over the number of PRIAM voxels
for rr = 1:size(vox_center,2)
    % Create a mask with all voxels that are inside the voxel
    mask = zeros(1,size(XYZ,2));
    sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
    distance2voxctr = sqrt(sum((XYZ-repmat([vox_center(1,rr) vox_center(2,rr) vox_center(3,rr)].',[1 size(XYZ,2)])).^2,1)); % LR-AP-FH
    sphere_mask(distance2voxctr <= sphere_radius) = 1;
    
    mask(sphere_mask == 1) = 1;
    XYZ_sphere = XYZ(:,sphere_mask == 1);
    
    tri = delaunayn([vox_corner(:,:,rr).'; [vox_center(1,rr) vox_center(2,rr) vox_center(3,rr)]]);
    tn = tsearchn([vox_corner(:,:,rr).'; [vox_center(1,rr) vox_center(2,rr) vox_center(3,rr)]], tri, XYZ_sphere.');
    isinside = ~isnan(tn);
    mask(sphere_mask==1) = isinside;
    
    % Take over the voxel dimensions from the structural
    mask = reshape(mask, vol_sos.dim);
    
    % Get the input file name
    [path,filename,~]   = fileparts(MRSCont.files{kk});
    % For batch analysis, get the last two sub-folders (e.g. site and
    % subject)
    path_split          = regexp(path,filesep,'split');
    if length(path_split) > 2
        saveName = [path_split{end-1} '_' path_split{end} '_' filename];
    end
    % Set up saving location
    

    saveDestination = fullfile(MRSCont.outputFolder,'SenseReconstruction', 'VoxelMasks');
    if ~exist(saveDestination,'dir')
        mkdir(saveDestination);
    end
    % Generate file name for the voxel mask NIfTI file to be saved under
    maskFile            = fullfile(saveDestination, [saveName '_Voxel' num2str(rr) '_VoxelMask.nii']);
        
    % Fill in the SPM volume header information
    vol_mask.fname   = maskFile;
    vol_mask.dim     = vol_sos.dim;
    vol_mask.dt      = vol_sos.dt;
    vol_mask.mat     = vol_sos.mat;
    vol_mask.descrip = 'MRS_voxel_mask';
    
    % Write the SPM volume to disk
    vol_mask = spm_write_vol(vol_mask,mask);
    
    % Store filename of NIfTI 
    MRSCont.SENSE{kk}.vol_mask{rr}  = vol_mask.fname;
end


% Store voxel center and corners for output figure
MRSCont.SENSE{kk}.sos_image = niiSOSFile;
MRSCont.SENSE{kk}.sense_image = niiSENSEFile;
MRSCont.SENSE{kk}.SOS_max   = SOS_max;
MRSCont.SENSE{kk}.voxel_ctr = vox_center;
MRSCont.SENSE{kk}.vox_corner = vox_corner;

% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

% Optional: 
if launchViewer
    nii_viewer(niiSOSFile, {MRSCont.SENSE{kk}.vol_mask{1}, MRSCont.SENSE{kk}.vol_mask{2}});
    nii_viewer(niiSENSEFile, {MRSCont.SENSE{kk}.vol_mask{1}, MRSCont.SENSE{kk}.vol_mask{2}});
end


% Generate SENSE matrix
niiSOS      = nii_tool('load', niiSOSFile);
niiSENSE    = nii_tool('load', niiSENSEFile);
% Loop over the number of PRIAM voxels
for rr = 1:size(vox_center,2)
    niiVox{rr} = nii_tool('load', MRSCont.SENSE{kk}.vol_mask{rr});
    for qq = 1:size(niiSENSE.img, 4)
        % Cut out the voxel mask from each coil image
        maskedCoilImage     = squeeze(niiSENSE.img(:,:,:,qq)) .* niiVox{rr}.img(:,:,:);
        % Calculate mean complex sensitivity across this voxel
        % NONZEROS ISN'T EXACTLY CORRECT BECAUSE YOU NULL THE ZEROS IN THE
        % VOXEL TOO - MAKE SURE YOU ACTUALLY CUT OUT THE MASK
        % IE DON'T MULTIPLY WITH THE MASK BUT SELECT THE INDICES THAT ARE 1
        % IN THE MASK
        %sensPerCoil(rr,qq)  = mean(nonzeros(maskedCoilImage));
        sensPerCoil(rr,qq)  = sum(sum(sum(maskedCoilImage))) ./ sum(sum(sum(niiVox{rr}.img)));
    end
end


% Get dimensions, orientation and global offset of the reference image scan
MRSCont.SENSE{kk}.refsc_sin_info   = get_sin_info([MRSCont.files_sense{kk}(1:end-4) '.sin']);
% Now pick only the coils that have been used in the actual MRS scan.
MRSCont.SENSE{kk}.nRecCoils = niiSENSE.hdr.dim(5); % save number of receiver coils
for ii = 1:MRSCont.SENSE{kk}.nRecCoils
    for ll = 1:length(MRSCont.SENSE{kk}.sin_info.channel_names)
        MRSCont.SENSE{kk}.CoilsUsedInMRS(ii) = strcmp(MRSCont.SENSE{kk}.refsc_sin_info.channel_names{ii},MRSCont.SENSE{kk}.sin_info.channel_names{ll});
        if MRSCont.SENSE{kk}.CoilsUsedInMRS(ii) == 1
            break
        end
    end
end

toDelete = find(MRSCont.SENSE{kk}.CoilsUsedInMRS==0);
sensPerCoil(:,toDelete) = [];
noise_array = MRSCont.SENSE{kk}.noise_array;
noise_array(toDelete,:) = [];



%%
% Normalize within each row of the coil combination coefficient matrix
% (separately for each voxel that we want to unfold)
for rr = 1:size(sensPerCoil,1)
    sensPerCoil(rr,:) = sensPerCoil(rr,:)./squeeze(sqrt(sum(abs(sensPerCoil(rr,:).^2))));
end

% Build sensitivity matrix S
S = sensPerCoil';

% Build noise correlation matrix
PSI = noise_array*noise_array';

% Build 
U = inv(S'*inv(PSI)*S)*S'*inv(PSI);
ga = inv(S'*inv(PSI)*S);
gb = S'*inv(PSI)*S;
MRSCont.SENSE{kk}.gfact = sqrt(ga.*gb);
MRSCont.SENSE{kk}.U = U;
MRSCont.SENSE{kk}.S = S;
MRSCont.SENSE{kk}.PSI = PSI;

disp(['gfactors: ' num2str(diag(real(MRSCont.SENSE{kk}.gfact))')])

end
