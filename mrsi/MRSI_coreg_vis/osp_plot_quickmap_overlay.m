function osp_plot_quickmap_overlay(quickmap, voxel_mask, vertices_display, nXvoxels, idx_start, idx_end, slices_per_row, tile_width, tile_height, custom_cmap, cmap_sz, alpha, orientation_info, display_info, isQC)
%% osp_plot_quickmap_overlay
%   Plots metabolite quickmap or metabolite map overlay on montage.
%
%   This function handles both continuous metabolite maps (quickmaps,
%   quantification results) and discrete QC maps with categorical coloring.
%
%   USAGE:
%       osp_plot_quickmap_overlay(quickmap, voxel_mask, vertices_display, nXvoxels, ...
%           idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
%           custom_cmap, cmap_sz, alpha, orientation_info, display_info, isQC)
%
%   ARGUMENTS:
%       quickmap        = 3D metabolite map data [nX, nY, nZ]
%       voxel_mask      = 3D binary mask for voxel selection [nX, nY, nZ]
%       vertices_display = Transformed vertices for display
%       nXvoxels        = Number of voxels in X dimension
%       idx_start       = Starting slice index
%       idx_end         = Ending slice index
%       slices_per_row  = Number of slices per row in montage
%       tile_width      = Width of each tile in montage
%       tile_height     = Height of each tile in montage
%       custom_cmap     = Colormap to use
%       cmap_sz         = Colormap scaling factor (max_val / 255)
%       alpha           = Transparency level (0-1)
%       orientation_info = Orientation information structure
%       display_info    = Display transformation information structure
%       isQC            = Flag for QC plots (0 = continuous, 1 = discrete/categorical)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2026-01-06: First version of the code.
%% Parse optional isQC parameter
if nargin < 15
    isQC = 0;
end

%% Transpose quickmap and voxel_mask to match vertex ordering
% Data is stored as [nX, nY, nZ] but vertex indexing assumes [nY, nX, nZ]
quickmap_t = permute(quickmap, [2, 1, 3]);
voxel_mask_t = permute(voxel_mask, [2, 1, 3]);

%% Check for single slice case
if ndims(voxel_mask_t) == 2
    isSingleSlice = 1;
    % Ensure 3D for consistent indexing
    voxel_mask_t = reshape(voxel_mask_t, [size(voxel_mask_t, 1), size(voxel_mask_t, 2), 1]);
    quickmap_t = reshape(quickmap_t, [size(quickmap_t, 1), size(quickmap_t, 2), 1]);
else
    isSingleSlice = 0;
end

%% Find voxels to plot
linear_indices = find(voxel_mask_t);

if isempty(linear_indices)
    warning('No voxels found in mask.');
    return;
end

%% Get voxel indices (now in [nY, nX, nZ] order after permute)
[yIdx, xIdx, zIdx] = ind2sub(size(voxel_mask_t), linear_indices);

%% Filter to slice range
valid_slices = (zIdx >= idx_start) & (zIdx <= idx_end);
xIdx = xIdx(valid_slices);
yIdx = yIdx(valid_slices);
zIdx = zIdx(valid_slices);

if isempty(xIdx)
    warning('No voxels in the specified slice range.');
    return;
end

%% Get colormap size for bounds checking
num_colors = size(custom_cmap, 1);

%% Plot each voxel
for vox = 1:length(xIdx)
    % Current voxel indices (1-indexed, in transposed space)
    x = xIdx(vox);
    y = yIdx(vox);
    z = zIdx(vox);
    
    % Get map value at this voxel (use transposed map)
    map_val = quickmap_t(y, x, z);
    
    % Handle color index based on QC mode
    if isQC
        % QC mode: discrete/categorical coloring
        % Values are expected to be integers (0, 1, 2, ...) representing categories
        if isnan(map_val)
            continue;
        end
        % Convert to 1-indexed colormap index
        map_ind = round(map_val) + 1;
        % Clamp to valid colormap range
        map_ind = max(1, min(num_colors, map_ind));
    else
        % Continuous mode: scale value to colormap
        % Skip if value is NaN, zero, or negative
        if isnan(map_val) || map_val <= 0
            continue;
        end
        % Convert to colormap index
        map_ind = round(map_val / cmap_sz) + 1;
        map_ind = max(1, min(num_colors-1, map_ind));
    end
    
    % Calculate vertex index (0-indexed coordinates)
    vertex_idx = ((y-1) * nXvoxels + (x-1)) * 8;
    
    % Calculate montage position based on slice
    slice_in_montage = z - idx_start;
    montage_col = mod(slice_in_montage, slices_per_row);
    montage_row = floor(slice_in_montage / slices_per_row);
    
    montage_x_offset = montage_col * tile_width;
    montage_y_offset = montage_row * tile_height;
    
    % Get voxel vertices (first 4 define the face)
    if vertex_idx >= 0 && vertex_idx + 4 <= size(vertices_display, 1)
        verts = vertices_display((1:4) + vertex_idx, :);
        
        % Consistent mapping: dim1 = Y (rows), dim2 = X (cols)
        plot_x = verts(:, 2) + montage_x_offset;
        plot_y = verts(:, 1) + montage_y_offset;
        
        % Create triangulation and plot
        try
            tri = delaunay(double(plot_x), double(plot_y));
            if isSingleSlice
                trisurf(tri, ...
                    plot_x, ...
                    plot_y, ...
                    zeros(4, 1), ...
                    'LineWidth', 1, ...
                    'FaceColor', custom_cmap(map_ind, :), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', alpha);
            else
                trisurf(tri, ...
                    plot_x, ...
                    plot_y, ...
                    verts(:, 3), ...
                    'LineWidth', 1, ...
                    'FaceColor', custom_cmap(map_ind, :), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', alpha);
            end
        catch
            % Skip voxels that fail triangulation (e.g., degenerate geometry)
            continue;
        end
    end
end

end