function [labels] = osp_get_orientation_labels(orientation_info, display_info)
%   This function determines anatomical orientation labels for each edge
%   of a displayed image after all spatial transformations have been applied.
%
%   The function accounts for permutations and flips applied during display
%   preparation to correctly identify what anatomical direction each edge
%   of the final displayed image represents.
%
%   Anatomical labels follow the RAS+ convention:
%       Positive directions: R (Right), A (Anterior), S (Superior)
%       Negative directions: L (Left), P (Posterior), I (Inferior)
%
%   USAGE:
%       labels = osp_get_orientation_labels(orientation_info, display_info);
%
%   INPUTS:
%       orientation_info = Struct from osp_analyze_orientation containing
%                          orientation analysis results.
%       display_info     = Struct from osp_prepare_display_slice containing:
%           .display_vert_world    - World axis for vertical display dimension
%           .display_horiz_world   - World axis for horizontal display dimension
%           .vert_sign_original    - Original sign of vertical mapping
%           .horiz_sign_original   - Original sign of horizontal mapping
%           .flipped_dims          - Array of dimensions that were flipped
%
%   OUTPUTS:
%       labels           = Struct containing orientation labels:
%           .left            - Anatomical label for left edge (e.g., 'R')
%           .right           - Anatomical label for right edge (e.g., 'L')
%           .top             - Anatomical label for top edge (e.g., 'A')
%           .bottom          - Anatomical label for bottom edge (e.g., 'P')
%           .slice_axis      - Slice axis name (e.g., 'S/I')
%           .slice_pos_label - Positive slice direction label (e.g., 'S')
%           .slice_neg_label - Negative slice direction label (e.g., 'I')
%           .vert_axis       - Vertical axis name (e.g., 'A/P')
%           .horiz_axis      - Horizontal axis name (e.g., 'L/R')
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    
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