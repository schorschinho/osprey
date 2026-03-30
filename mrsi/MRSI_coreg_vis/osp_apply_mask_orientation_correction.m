function mask_corrected = osp_apply_mask_orientation_correction(mask, orientation_info, display_info)
%   This function applies the same spatial transformations to masks that
%   were applied to the corresponding image data during display preparation.
%
%   This ensures that masks (e.g., voxel masks, segmentation masks) align
%   correctly with the displayed image after any permutation and flipping
%   operations that were performed for proper radiological display.
%
%   The function handles masks of different dimensionalities (2D or 3D)
%   and applies transformations only to valid dimensions.
%
%   USAGE:
%       mask_corrected = osp_apply_mask_orientation_correction(mask, orientation_info, display_info);
%
%   INPUTS:
%       mask             = 2D or 3D mask array to be transformed.
%       orientation_info = Struct from osp_analyze_orientation containing
%                          orientation analysis results.
%       display_info     = Struct from osp_prepare_display_slice containing:
%           .permuted       - Boolean indicating if permutation was applied
%           .permute_order  - [1x3] permutation order used
%           .flipped_dims   - Array of dimensions that were flipped
%
%   OUTPUTS:
%       mask_corrected   = Transformed mask aligned with displayed image.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    
    mask_corrected = mask;
    mask_size = size(mask);
    n_dims = length(mask_size);
    
    % First permute (same as image)
    if display_info.permuted
        perm = display_info.permute_order;
        
        if n_dims < 3
            perm_truncated = perm(perm <= n_dims);
            if length(perm_truncated) >= 2
                mask_corrected = permute(mask_corrected, perm_truncated);
            end
        else
            mask_corrected = permute(mask_corrected, perm);
        end
    end
    
    % Then flip (same dimensions as image)
    for flip_dim = display_info.flipped_dims
        if flip_dim <= ndims(mask_corrected)
            mask_corrected = flip(mask_corrected, flip_dim);
        end
    end
end