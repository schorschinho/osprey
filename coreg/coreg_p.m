% coreg_p.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [vol_mask] = coreg_p(in, dcm_folder, maskFile);
%
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions
% as the SPM volume containing a structural image. The voxel mask will be
% created with the geometry information stored in the FID-A structure "in".
% This routine works for the GE P data type.

% INPUTS:
% in        = Input data structure.
% dcm_folder = Path to the folder containing DICOM images of the structural
%           image that the voxel defined in the 'geometry' field of the
%           input data structure 'in' is supposed to be co-registered to.
% maskFile  = Filename under which the SPM volume of the co-registered
%           voxel mask is supposed to be saved.
%
% OUTPUTS:
% vol_mask  = SPM volume of the coregistered voxel mask.
% T1_max    = maximum intensity of the image volume.
% vol_image	= SPM volume of the structural image after DICOM conversion.

function [vol_mask, T1_max, vol_image, voxel_ctr] = coreg_p(in, dcm_folder, maskFile)

% Deactivate MATLAB warnings and load geometry parameters
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');
geom = in.geometry;

%%% 1. PREPARE THE MRS VOXEL COORDINATES %%%
% Convert from RAS to LPS
VoxOffs  = geom.pos .* [-1 -1 1];
tlhc_LPS = geom.rot.tlhc .* [-1 -1 1];
trhc_LPS = geom.rot.trhc .* [-1 -1 1];
brhc_LPS = geom.rot.brhc .* [-1 -1 1];

e1_SVS_n = trhc_LPS - tlhc_LPS;
e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
e2_SVS_n = brhc_LPS - trhc_LPS;
e2_SVS_n = e2_SVS_n ./ norm(e2_SVS_n);
e3_SVS_n = -cross(e1_SVS_n, e2_SVS_n);

[~,orientation_SVS] = max(abs(e3_SVS_n));

if orientation_SVS == 3     % axial
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e2_SVS_n;
    e3_SVS_n2 = e3_SVS_n;
elseif orientation_SVS == 2 % coronal
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e3_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
elseif orientation_SVS == 1 % sagittal
    e1_SVS_n2 = e3_SVS_n;
    e2_SVS_n2 = e1_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
end

angulations = get_euler(e1_SVS_n2, e2_SVS_n2, e3_SVS_n2);

e1_SVS = geom.size.dim1 * e1_SVS_n2;
e2_SVS = geom.size.dim2 * e2_SVS_n2;
e3_SVS = geom.size.dim3 * e3_SVS_n2;

% LPS gives center of voxel
LPS_SVS_edge = VoxOffs - 0.5 * e1_SVS ...
                        - 0.5 * e2_SVS ...
                        - 0.5 * e3_SVS;


%%% 2. PREPARE THE STRUCTURAL IMAGE
% Read all DICOM files into one volume
dcm_list = dir(dcm_folder);
dcm_list = dcm_list(~ismember({dcm_list.name}, {'.','..','.DS_Store'}));
hidden = logical(ones(1,length(dcm_list)));
for jj = 1:length(dcm_list) 
    if strcmp(dcm_list(jj).name(1),'.')
        hidden(jj) = 0;
    end
end
dcm_list = dcm_list(hidden);%delete hidden files 
dcm_list = cellstr(char(dcm_list.name));
dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.nii'))); %#ok<*STRCLFH>
dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.mat')));
for jj = 1:length(dcm_list)
    dcm_list{jj} = [dcm_folder filesep dcm_list{jj}];
end
dcm_hdr = spm_dicom_headers(char(dcm_list));
nii_file_dir = spm_dicom_convert_osp(dcm_hdr, 'all', 'flat', 'nii', fullfile(dcm_folder)); % create NIFTI file of T1 image
nii_file = nii_file_dir.files{1};
% Create SPM volume and read in the NIfTI file with the structural image.
vol_image = spm_vol(nii_file);
[T1,~]  = spm_read_vols(vol_image);
T1_max    = max(T1(:));

slice_location = zeros(1,length(dcm_list));
for jj = 1:length(dcm_list)
    slice_location(jj) = dcm_hdr{jj}.SliceLocation;
end

% Order slices according to slice position
[~,order_index] = sort(slice_location);
tmp = dcm_hdr;
for jj = 1:length(dcm_list)
    dcm_hdr{jj} = tmp{order_index(jj)};
end

MRI_voxel_size = [dcm_hdr{1}.PixelSpacing(1) ...
                  dcm_hdr{1}.PixelSpacing(2) ...
                  dcm_hdr{1}.SpacingBetweenSlices];

MRI_dim = [dcm_hdr{1}.Rows ...
           dcm_hdr{1}.Columns ...
           dcm_hdr{1}.ImagesInAcquisition];

e1_MRI_n = dcm_hdr{1}.ImageOrientationPatient(1:3);
e2_MRI_n = dcm_hdr{1}.ImageOrientationPatient(4:6);
e3_MRI_n = cross(e1_MRI_n, e2_MRI_n); % e3 vector is perpendicular to the slice orientation

[~,orientation_MRI] = max(abs(e3_MRI_n));
if orientation_MRI == 2 % coronal
    e3_MRI_n = -e3_MRI_n;
end

if orientation_MRI == 3     % axial
    e1_MRI_n2 = e1_MRI_n;
    e2_MRI_n2 = e2_MRI_n;
    e3_MRI_n2 = e3_MRI_n;
    MRI_voxel_size = MRI_voxel_size([2 1 3]);
    MRI_dim = MRI_dim([2 1 3]);
elseif orientation_MRI == 2 % coronal
    e1_MRI_n2 = e1_MRI_n;
    e2_MRI_n2 = e3_MRI_n;
    e3_MRI_n2 = e2_MRI_n;
    MRI_voxel_size = MRI_voxel_size([2 3 1]);
    MRI_dim = MRI_dim([2 3 1]);
elseif orientation_MRI == 1 % sagittal
    e1_MRI_n2 = e3_MRI_n;
    e2_MRI_n2 = e1_MRI_n;
    e3_MRI_n2 = e2_MRI_n;
    MRI_voxel_size = MRI_voxel_size([3 2 1]);
    MRI_dim = MRI_dim([3 2 1]);
end

% LPS_edge gives location of the edge of the image volume
LPS_MRI_center = dcm_hdr{1}.ImagePositionPatient;
LPS_MRI_edge = LPS_MRI_center - 0.5 * MRI_voxel_size(1) * e1_MRI_n2 ...
                              - 0.5 * MRI_voxel_size(2) * e2_MRI_n2 ...
                              - 0.5 * MRI_voxel_size(3) * e3_MRI_n2;

% Create voxel mask
E_MRI = [e1_MRI_n2 e2_MRI_n2 e3_MRI_n2];
c_MRS = VoxOffs';
c_MRI = E_MRI' * (c_MRS - LPS_MRI_edge);
d_MRI = c_MRI ./ MRI_voxel_size';
s_MRS = sqrt(sum([geom.size.dim1 .^2, geom.size.dim2 .^2, geom.size.dim3 .^2]))/2;
d_MRS = s_MRS ./ MRI_voxel_size';

[Xm,Ym,Zm] = ndgrid(1:MRI_dim(1), 1:MRI_dim(2), 1:MRI_dim(3));
X = LPS_MRI_center(1) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(1) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(1) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(1);
Y = LPS_MRI_center(2) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(2) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(2) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(2);
Z = LPS_MRI_center(3) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(3) + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(3) + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(3);

P_1 = LPS_SVS_edge;
P_2 = LPS_SVS_edge + e1_SVS; % L
P_3 = LPS_SVS_edge + e2_SVS; % P
P_4 = LPS_SVS_edge + e3_SVS; % S
A = zeros(3,1);
mask = zeros(MRI_dim);

for e1 = max(floor(d_MRI(1) - d_MRS(1)), 0) : min(ceil(d_MRI(1) + d_MRS(1)), size(mask,1))         % L
    for e2 = max(floor(d_MRI(2) - d_MRS(2)), 0) : min(ceil(d_MRI(2) + d_MRS(2)), size(mask,2))     % P
        for e3 = max(floor(d_MRI(3) - d_MRS(3)), 0) : min(ceil(d_MRI(3) + d_MRS(3)), size(mask,3)) % S
            A(1) = X(e1,e2,e3);
            A(2) = Y(e1,e2,e3);
            A(3) = Z(e1,e2,e3);
            % Distance of A to planes in SI direction
            d_5 = e3_SVS_n2 * (A - P_1');
            d_6 = -e3_SVS_n2 * (A - P_4');
            if d_5 >= 0 && d_6 >= 0
                % Distance of A to planes in AP direction
                d_3 = e2_SVS_n2 * (A - P_1');
                d_4 = -e2_SVS_n2 * (A - P_3');
                if d_3 >= 0 && d_4 >= 0
                    % Distance of A to planes in RL direction
                    d_1 = e1_SVS_n2 * (A - P_1');
                    d_2 = -e1_SVS_n2 * (A - P_2');
                    if d_1 >= 0 && d_2 >= 0
                        mask(e1,e2,e3) = 1;
                    end
                end
            end
        end
    end
end

if orientation_MRI == 2     % coronal
    mask = permute(mask, [3 1 2]);
    mask = flip(mask,3);
elseif orientation_MRI == 1 % sagittal
    mask = permute(mask, [2 3 1]);
end
mask = flip(mask,2);

% Fill in the SPM volume header information
vol_mask.fname   = maskFile;
vol_mask.dim     = vol_image.dim;
vol_mask.dt      = vol_image.dt;
vol_mask.mat     = vol_image.mat;
vol_mask.descrip = 'MRS_voxel_mask';

% Write the SPM volume to disk
vol_mask = spm_write_vol(vol_mask,mask);

% Store voxel centre for output figure
VoxOffs(1:2) = -VoxOffs(1:2);
voxel_ctr = VoxOffs;

% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end



function euler_angles = get_euler(r1, r2, r3)

r1(3) = -r1(3);
r2(3) = -r2(3);
r3(3) = -r3(3);

if abs(r3(1)) ~= 1
    theta1 = -asin(r3(1));
    %theta2 = pi - theta1;
    psi1 = atan2(r3(2)/cos(theta1), r3(3)/cos(theta1));
    %psi2 = atan2(r3(2)/cos(theta2), r3(3)/cos(theta2));
    phi1 = atan2(r2(1)/cos(theta1), r1(1)/cos(theta1));
    %phi2 = atan2(r2(1)/cos(theta2), r1(1)/cos(theta2));
else
    phi1 = 0;
    if r3(1) == -1
        theta1 = pi/2;
        psi1 = phi1 + atan2(r1(2), r1(3));
    else
        theta1 = -pi/2;
        psi1 = -phi1 + atan2(-r1(2), -r1(3));
    end
end

euler_angles(1) = round(-phi1*180/pi);
euler_angles(2) = round(-psi1*180/pi);
euler_angles(3) = round(theta1*180/pi);

end
