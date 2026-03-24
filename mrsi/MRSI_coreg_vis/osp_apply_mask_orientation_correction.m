function mask_corrected = osp_apply_mask_orientation_correction(mask, orientation_info, display_info)
    % Apply the same transformations to masks as were applied to the image
    %
    % This ensures masks align correctly with the displayed image
    
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