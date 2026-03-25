function resize_img(imnames, ref_img, ismask,interpolate)
% resize_img -- Resample images to match the affine matrix, resolution, and bounding box of a reference image
% resize_img(imnames, ref_img, ismask)
%
% - imnames: List of input image filenames to be resized.
% - ref_img: Reference image whose affine, resolution, and bounding box will be used.
% - ismask: Boolean, whether input images are binary masks (avoid interpolation effects on masks).
%
% The output images will be prefixed with 'r', and they will match the
% voxel dimensions, affine matrix, and bounding box of the reference image.

if nargin < 3
    ismask = false;
end
if nargin < 4
    interpolate = [];
end

% Load reference image
refVol = spm_vol(ref_img);
refDim = refVol.dim;          % Dimensions of the reference image
if ~isempty(interpolate)
    refDim(1) = refDim(1)*interpolate;
    refDim(2) = refDim(2)*interpolate;
end
refMat = refVol.mat;          % Affine matrix of the reference image
if ~isempty(interpolate)
    scalingFactors = [interpolate interpolate 1];
    refMat(1:3,1:3) = refMat(1:3,1:3) * diag(1 ./ scalingFactors);   
end
refRes = abs(diag(refMat(1:3, 1:3)))';  % Voxel resolution from affine matrix

% Resample each target image
vols = spm_vol(imnames);
for V = vols'
    % Update output filename
    [pth, nam, ext] = fileparts(V.fname);
    rFilename = fullfile(pth, ['r', nam, ext]);
    
    % Create new volume for the output image
    VO = V;
    VO.fname = rFilename;
    VO.dim = refDim;
    VO.mat = refMat;

    % Create new image in memory
    VO = spm_create_vol(VO);
    spm_progress_bar('Init', refDim(3), 'Reslicing...', 'Planes completed');

    for i = 1:refDim(3)
        % Transformation matrix to align input image to reference space
        M = inv(spm_matrix([0 0 -i]) * inv(VO.mat) * V.mat);
        
        % Reslice image with linear interpolation (or nearest-neighbor for masks)
        interp = 0; % Nearest neighbor for masks
        if ~ismask
            interp = 1; % Linear interpolation for non-mask images
        end
        img = spm_slice_vol(V, M, refDim(1:2), interp);
        
        % Round mask values to avoid interpolation artifacts
        if ismask
            img = round(img);
        end
        
        % Write the image slice to the output file
        spm_write_plane(VO, img, i);
        spm_progress_bar('Set', i);
    end
    spm_progress_bar('Clear');
end
end