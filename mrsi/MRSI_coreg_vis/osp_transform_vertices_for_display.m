function [vertices_display] = osp_transform_vertices_for_display(vertices_voxel, img_original, orientation_info, display_info, MRSI_vol)
    % Transform vertex coordinates to match the displayed image
    %
    % Applies the same transformations to vertices as were applied to the image:
    %   1. Permutation (same order as image)
    %   2. Flips (same dimensions as image)
    %   3. Half-voxel offset (for non-axial orientations)
    %
    % After transformation:
    %   vertices_display(:,1) = row index in displayed image (plot Y)
    %   vertices_display(:,2) = column index in displayed image (plot X)
    %   vertices_display(:,3) = slice index
    
    img_size = size(img_original);
    
    % === Calculate scale factors for half-voxel offset ===
    if nargin >= 5 && ~isempty(MRSI_vol)
        mrsi_dim = MRSI_vol.dim(1:2);
        display_dim = img_size(1:2);
        scale_factor = display_dim ./ mrsi_dim;
        half_voxel_offset = scale_factor * 0.5;
    else
        half_voxel_offset = [0.5, 0.5];
    end
    
    % === Apply permutation (same as image) ===
    if display_info.permuted
        perm = display_info.permute_order;
        n_dims = length(img_size);
        
        % Extend half_voxel_offset to 3D for permutation
        half_voxel_offset_3d = [half_voxel_offset(:)', 0];
        
        if n_dims < 3
            perm_truncated = perm(perm <= n_dims);
            vertices_display = vertices_voxel(:, perm_truncated);
            img_size = img_size(perm(1:n_dims));
            half_voxel_offset = half_voxel_offset_3d(perm(1:n_dims));
        else
            vertices_display = vertices_voxel(:, perm);
            img_size = img_size(perm);
            half_voxel_offset_3d = half_voxel_offset_3d(perm);
            half_voxel_offset = half_voxel_offset_3d(1:2);
        end
    else
        vertices_display = vertices_voxel;
    end
    
    % === Apply flips (same dimensions as image) ===
    for flip_dim = display_info.flipped_dims
        if flip_dim <= length(img_size)
            vertices_display(:, flip_dim) = img_size(flip_dim) - vertices_display(:, flip_dim) + 1;
        end
    end
    
    % === Apply half-voxel offset ===
    % Needed for non-axial orientations where coordinate systems don't align perfectly
    switch orientation_info.slice_orientation
        case 'axial'
            % No offset needed for axial
        case {'sagittal', 'coronal', 'oblique'}
            if display_info.permuted
                vertices_display(:, 1) = vertices_display(:, 1);
                vertices_display(:, 2) = vertices_display(:, 2);
            end
    end
end
