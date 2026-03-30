function [orientation_info] = osp_analyze_orientation(AffineMat)
%   This function analyzes a NIfTI affine transformation matrix to determine
%   the image orientation relative to the RAS+ world coordinate system.
%
%   World coordinate system (RAS+):
%       X (1): Left (-) to Right (+)
%       Y (2): Posterior (-) to Anterior (+)
%       Z (3): Inferior (-) to Superior (+)
%
%   The function determines slice orientation (axial, sagittal, coronal, or
%   oblique) and provides mappings between image dimensions and world axes
%   for proper display and spatial referencing.
%
%   USAGE:
%       orientation_info = osp_analyze_orientation(AffineMat);
%
%   INPUTS:
%       AffineMat   = 4x4 affine transformation matrix from NIfTI header.
%
%   OUTPUTS:
%       orientation_info = Struct containing:
%           .AffineMat           - Original affine matrix
%           .dir_cos             - 3x3 direction cosines matrix
%           .vox_size            - [1x3] voxel dimensions
%           .img_dim_to_world    - [1x3] world axis for each image dimension
%           .img_dim_alignment   - [1x3] alignment strength (0-1)
%           .world_to_img_dim    - [1x3] image dimension for each world axis
%           .dim_signs           - [1x3] sign of each dimension mapping
%           .slice_orientation   - 'axial', 'sagittal', 'coronal', or 'oblique'
%           .standard_horiz_world - World axis for horizontal display
%           .standard_vert_world  - World axis for vertical display
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on NIfTI orientation conventions as described in
%       the NIfTI-1 data format specification.
%       https://nifti.nimh.nih.gov/nifti-1
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    
    R = AffineMat(1:3, 1:3);
    vox_size = sqrt(sum(R.^2, 1));
    dir_cos = R ./ vox_size;
    
    orientation_info.AffineMat = AffineMat;
    orientation_info.dir_cos = dir_cos;
    orientation_info.vox_size = vox_size;
    
    % For each image dimension, find dominant world axis
    [~, img_dim_to_world] = max(abs(dir_cos), [], 1);
    orientation_info.img_dim_to_world = img_dim_to_world;
    
    % For each world axis, find corresponding image dimension
    world_to_img_dim = zeros(1, 3);
    for w = 1:3
        [~, world_to_img_dim(w)] = max(abs(dir_cos(w, :)));
    end
    orientation_info.world_to_img_dim = world_to_img_dim;
    
    % Get signs for each image dimension's mapping to world axis
    dim_signs = zeros(1, 3);
    for d = 1:3
        dim_signs(d) = sign(dir_cos(img_dim_to_world(d), d));
    end
    orientation_info.dim_signs = dim_signs;
    
    % Determine slice orientation based on 3rd image dimension (slice direction)
    slice_world_axis = img_dim_to_world(3);
    oblique_threshold = 0.9;
    
    max_alignment = max(abs(dir_cos(:, 3)));
    
    if max_alignment < oblique_threshold
        orientation_info.slice_orientation = 'oblique';
    else
        switch slice_world_axis
            case 1  % Slices perpendicular to L-R
                orientation_info.slice_orientation = 'sagittal';
            case 2  % Slices perpendicular to A-P
                orientation_info.slice_orientation = 'coronal';
            case 3  % Slices perpendicular to S-I
                orientation_info.slice_orientation = 'axial';
        end
    end
    
    % Determine standard display arrangement for each orientation
    % Standard views (for display):
    %   Axial: horizontal=L-R (1), vertical=A-P (2)
    %   Sagittal: horizontal=A-P (2), vertical=S-I (3)
    %   Coronal: horizontal=L-R (1), vertical=S-I (3)
    
    switch orientation_info.slice_orientation
        case 'axial'
            orientation_info.standard_horiz_world = 1;  % L-R
            orientation_info.standard_vert_world = 2;   % A-P
        case 'sagittal'
            orientation_info.standard_horiz_world = 2;  % A-P
            orientation_info.standard_vert_world = 3;   % S-I
        case 'coronal'
            orientation_info.standard_horiz_world = 1;  % L-R
            orientation_info.standard_vert_world = 3;   % S-I
        case 'oblique'
            % Use the first two in-plane dimensions
            orientation_info.standard_horiz_world = img_dim_to_world(2);
            orientation_info.standard_vert_world = img_dim_to_world(1);
    end
    
    % fprintf('Detected orientation: %s\n', orientation_info.slice_orientation);
    % fprintf('Image dim to world: [%d, %d, %d]\n', img_dim_to_world);
    % fprintf('World to image dim: [%d, %d, %d]\n', world_to_img_dim);
end