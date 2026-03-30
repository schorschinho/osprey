function [img_display, display_info] = osp_prepare_image_for_display(img, orientation_info, convention)
%% [img_display, display_info] = osp_prepare_image_for_display(img, orientation_info, convention)
%   This function prepares a 3D image volume for display with correct
%   anatomical orientation by applying appropriate permutations and flips.
%
%   The function ensures the image is oriented according to standard
%   radiological or neurological display conventions:
%       1. Permutes dimensions so that:
%          - Dimension 1 (rows) = vertical display axis
%          - Dimension 2 (cols) = horizontal display axis
%          - Dimension 3 = slices
%       2. Flips dimensions for proper anatomical orientation
%
%   Standard display conventions by slice orientation:
%       Axial:    Horizontal=L-R, Vertical=A-P, Anterior at top
%       Sagittal: Horizontal=A-P, Vertical=S-I, Anterior at left, Superior at top
%       Coronal:  Horizontal=L-R, Vertical=S-I, Superior at top
%
%
%   USAGE:
%       [img_display, display_info] = osp_prepare_image_for_display(img, orientation_info);
%       [img_display, display_info] = osp_prepare_image_for_display(img, orientation_info, convention);
%
%   INPUTS:
%       img              = 2D or 3D image volume array.
%       orientation_info = Struct from osp_analyze_orientation containing
%                          orientation analysis results.
%       convention       = Display convention string (optional):
%                          'radiological' (default) or 'neurological'
%
%   OUTPUTS:
%       img_display      = Transformed image ready for display.
%       display_info     = Struct containing transformation details:
%           .convention            - Display convention used
%           .img_size_original     - Original image dimensions
%           .img_size_display      - Final display image dimensions
%           .flipped_dims          - Array of dimensions that were flipped
%           .permute_order         - [1x3] permutation order applied
%           .permuted              - Boolean indicating if permutation applied
%           .is_single_slice       - Boolean indicating single slice input
%           .display_vert_world    - World axis for vertical display dimension
%           .display_horiz_world   - World axis for horizontal display dimension
%           .vert_sign_original    - Original sign of vertical axis mapping
%           .horiz_sign_original   - Original sign of horizontal axis mapping
%           .vert_img_dim_original - Original image dim that became rows
%           .horiz_img_dim_original - Original image dim that became columns
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    
    if nargin < 3
        convention = 'radiological';
    end
    
    display_info.convention = convention;
    display_info.img_size_original = size(img);
    display_info.flipped_dims = [];
    display_info.permute_order = [1, 2, 3];
    display_info.permuted = false;
    display_info.is_single_slice = false;
    
    % Check if single slice
    img_size = size(img);
    n_dims = length(img_size);
    
    if n_dims < 3 || img_size(3) == 1
        display_info.is_single_slice = true;
        
        % Ensure 3D for consistent processing
        if n_dims < 3
            img = reshape(img, [img_size, 1]);
            img_size = size(img);
        end
    end
    
    dir_cos = orientation_info.dir_cos;
    world_to_img_dim = orientation_info.world_to_img_dim;
    
    % Get the standard arrangement for this orientation
    std_horiz = orientation_info.standard_horiz_world;
    std_vert = orientation_info.standard_vert_world;
    
    % Find which image dimensions currently correspond to these world axes
    current_horiz_dim = world_to_img_dim(std_horiz);
    current_vert_dim = world_to_img_dim(std_vert);
    slice_dim = 3;
    
    % Build permutation order: [new_dim1, new_dim2, new_dim3]
    % new_dim1 = rows (vertical in display)
    % new_dim2 = columns (horizontal in display)
    % new_dim3 = slices
    permute_order = [current_vert_dim, current_horiz_dim, slice_dim];
    
    % Check if permutation is needed and valid
    needs_permutation = ~isequal(permute_order, [1, 2, 3]);
    is_valid_permutation = length(unique(permute_order)) == 3;
    
    if needs_permutation && is_valid_permutation
        display_info.permuted = true;
        display_info.permute_order = permute_order;
        img_display = permute(img, permute_order);
        img_size = size(img_display);
        % fprintf('Permuting image: [%d, %d, %d]\n', permute_order);
    else
        img_display = img;
        if ~is_valid_permutation && needs_permutation
            warning('Invalid permutation order: [%d, %d, %d]. Using original.', permute_order);
        end
    end
    
    % Store the world axes for the display dimensions
    display_info.display_vert_world = std_vert;
    display_info.display_horiz_world = std_horiz;
    
    % Get the direction signs for the display dimensions
    % After permutation, dim1 = vertical, dim2 = horizontal
    if display_info.permuted
        perm = display_info.permute_order;
        vert_img_dim_original = perm(1);  % Original dim that became rows
        horiz_img_dim_original = perm(2); % Original dim that became cols
    else
        vert_img_dim_original = 1;
        horiz_img_dim_original = 2;
    end
    
    vert_sign = sign(dir_cos(std_vert, vert_img_dim_original));
    horiz_sign = sign(dir_cos(std_horiz, horiz_img_dim_original));
    
    display_info.vert_sign_original = vert_sign;
    display_info.horiz_sign_original = horiz_sign;
    display_info.vert_img_dim_original = vert_img_dim_original;
    display_info.horiz_img_dim_original = horiz_img_dim_original;
    
    % Apply flips based on orientation and convention
    switch orientation_info.slice_orientation
        case 'axial'
            % Vertical = A-P: A should be at top (low row index)
            % If vert_sign > 0: increasing index = toward A, A at high index, need flip
            if vert_sign > 0
                img_display = flip(img_display, 1);
                display_info.flipped_dims = [display_info.flipped_dims, 1];
            end
            
            % Horizontal = L-R: depends on convention
            if strcmp(convention, 'radiological')
                % Radiological: L on right (high column index)
                % If horiz_sign > 0: R at high index, need flip to put L at high index
                if horiz_sign > 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            else  % neurological
                % Neurological: L on left (low column index)
                % If horiz_sign < 0: L at high index, need flip to put L at low index
                if horiz_sign < 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            end
            
        case 'sagittal'
            % Vertical = S-I: S should be at top (low row index)
            % If vert_sign > 0: S at high index, need flip
            if vert_sign > 0
                img_display = flip(img_display, 1);
                display_info.flipped_dims = [display_info.flipped_dims, 1];
            end
            
            % Horizontal = A-P: A should be at left (low column index)
            % If horiz_sign > 0: A at high index, need flip
            if horiz_sign > 0
                img_display = flip(img_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
            
        case 'coronal'
            % Vertical = S-I: S should be at top
            if vert_sign > 0
                img_display = flip(img_display, 1);
                display_info.flipped_dims = [display_info.flipped_dims, 1];
            end
            
            % Horizontal = L-R: depends on convention
            if strcmp(convention, 'radiological')
                if horiz_sign > 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            else
                if horiz_sign < 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            end
            
        case 'oblique'
            % Best effort for oblique
            si_in_vert = dir_cos(3, vert_img_dim_original);
            if si_in_vert > 0
                img_display = flip(img_display, 1);
                display_info.flipped_dims = [display_info.flipped_dims, 1];
            end
            
            lr_in_horiz = dir_cos(1, horiz_img_dim_original);
            if strcmp(convention, 'radiological')
                if lr_in_horiz > 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            else
                if lr_in_horiz < 0
                    img_display = flip(img_display, 2);
                    display_info.flipped_dims = [display_info.flipped_dims, 2];
                end
            end
    end
    
    display_info.img_size_display = size(img_display);
    
    % fprintf('Permuted: %d, Flipped dims: [%s]\n', display_info.permuted, num2str(display_info.flipped_dims));
    % fprintf('Display image size: [%d, %d, %d]\n', size(img_display, 1), size(img_display, 2), size(img_display, 3));
end