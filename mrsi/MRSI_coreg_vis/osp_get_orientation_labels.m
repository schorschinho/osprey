function [labels] = osp_get_orientation_labels(orientation_info, display_info)
    % Determine orientation labels for each edge of the displayed image
    %
    % After all transformations (permute + flip), determine what anatomical
    % direction each edge of the image represents
    %
    % Returns struct with:
    %   labels.left, labels.right: horizontal axis labels
    %   labels.top, labels.bottom: vertical axis labels
    %   labels.slice_axis: which world axis slices are along
    
    pos_labels = {'R', 'A', 'S'};  % Positive direction: Right, Anterior, Superior
    neg_labels = {'L', 'P', 'I'};  % Negative direction: Left, Posterior, Inferior
    axis_names = {'L/R', 'A/P', 'S/I'};
    
    % Get the world axes for display dimensions
    vert_world = display_info.display_vert_world;
    horiz_world = display_info.display_horiz_world;
    
    % Determine slice world axis
    slice_world_axis = orientation_info.img_dim_to_world(3);
    
    % Determine the effective sign after all transformations
    vert_sign = display_info.vert_sign_original;
    horiz_sign = display_info.horiz_sign_original;
    
    % Account for flips
    if ismember(1, display_info.flipped_dims)
        vert_sign = -vert_sign;
    end
    if ismember(2, display_info.flipped_dims)
        horiz_sign = -horiz_sign;
    end
    
    % Vertical axis (rows): row 1 at top, row N at bottom
    % If vert_sign > 0 after transforms: increasing row = positive world direction
    if vert_sign > 0
        labels.top = neg_labels{vert_world};
        labels.bottom = pos_labels{vert_world};
    else
        labels.top = pos_labels{vert_world};
        labels.bottom = neg_labels{vert_world};
    end
    
    % Horizontal axis (cols): col 1 at left, col N at right
    % If horiz_sign > 0 after transforms: increasing col = positive world direction
    if horiz_sign > 0
        labels.left = neg_labels{horiz_world};
        labels.right = pos_labels{horiz_world};
    else
        labels.left = pos_labels{horiz_world};
        labels.right = neg_labels{horiz_world};
    end
    
    % Slice axis information
    labels.slice_axis = axis_names{slice_world_axis};
    labels.slice_pos_label = pos_labels{slice_world_axis};
    labels.slice_neg_label = neg_labels{slice_world_axis};
    
    labels.vert_axis = axis_names{vert_world};
    labels.horiz_axis = axis_names{horiz_world};
end