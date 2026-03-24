function [orientation_info] = osp_analyze_orientation(AffineMat)
    % Analyze the affine matrix to determine image orientation
    %
    % World coordinate system (RAS+):
    %   X: Left (-) to Right (+)
    %   Y: Posterior (-) to Anterior (+)
    %   Z: Inferior (-) to Superior (+)
    %
    % Returns orientation_info struct with:
    %   - slice_orientation: 'axial', 'sagittal', 'coronal', or 'oblique'
    %   - dir_cos: direction cosines matrix
    %   - vox_size: voxel dimensions
    %   - img_dim_to_world: which world axis each image dim maps to
    %   - world_to_img_dim: which image dim each world axis maps to
    %   - standard_horiz_world: which world axis should be horizontal
    %   - standard_vert_world: which world axis should be vertical
    
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