function osp_plot_tissue_overlay(tissue_map, brain_mask, vertices_display, nXvoxels, idx_start, idx_end, slices_per_row, tile_width, tile_height, custom_cmap, cmap_sz, alpha, orientation_info, display_info)
%% osp_plot_tissue_overlay
%   Plots tissue segmentation map overlay on montage.
%
%   The tissue maps from MRSCont.seg.tissue are in MRSI voxel space [nX, nY, nZ].
%   The vertex indexing uses: vertex_idx = (y * nXvoxels + x) * 8
%   So we need to swap dimensions 1 and 2 to match.
%
%   USAGE:
%       osp_plot_tissue_overlay(tissue_map, brain_mask, vertices_display, nXvoxels, ...
%           idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
%           custom_cmap, cmap_sz, alpha, orientation_info, display_info)
%
%   ARGUMENTS:
%       tissue_map      = 3D tissue map in MRSI voxel space [nX, nY, nZ]
%       brain_mask      = 3D brain mask in MRSI voxel space [nX, nY, nZ]
%       vertices_display = transformed vertex coordinates for display
%       nXvoxels        = number of X voxels in MRSI grid
%       idx_start       = start slice index
%       idx_end         = end slice index
%       slices_per_row  = number of slices per row in montage
%       tile_width      = width of each tile in pixels
%       tile_height     = height of each tile in pixels
%       custom_cmap     = colormap to use (256 colors)
%       cmap_sz         = colormap step size (typically 1/255)
%       alpha           = transparency level (0-1)
%       orientation_info = orientation analysis struct
%       display_info    = display preparation struct
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2026-01-06: First version of the code.
%% Transpose tissue map and brain mask to match vertex ordering
% Tissue maps are stored as [nX, nY, nZ] but vertex indexing assumes [nY, nX, nZ]
% vertex_idx = (y * nXvoxels + x) * 8, where y is outer loop, x is inner
tissue_map_t = permute(tissue_map, [2, 1, 3]);  % Now [nY, nX, nZ]
brain_mask_t = permute(brain_mask, [2, 1, 3]);  % Now [nY, nX, nZ]

%% Find brain voxels to plot
linear_indices = find(brain_mask_t);

if isempty(linear_indices)
    warning('No brain voxels found in mask.');
    return;
end

%% Get voxel indices (now in [nY, nX, nZ] order after permute)
[yIdx, xIdx, zIdx] = ind2sub(size(brain_mask_t), linear_indices);

%% Filter to slice range
valid_slices = (zIdx >= idx_start) & (zIdx <= idx_end);
xIdx = xIdx(valid_slices);
yIdx = yIdx(valid_slices);
zIdx = zIdx(valid_slices);

if isempty(xIdx)
    warning('No voxels in the specified slice range.');
    return;
end

%% Plot each voxel
for vox = 1:length(xIdx)
    % Current voxel indices (1-indexed, in transposed space)
    x = xIdx(vox);  % X index (1 to nXvoxels)
    y = yIdx(vox);  % Y index (1 to nYvoxels)
    z = zIdx(vox);  % Z index (slice)
    
    % Get tissue map value at this voxel (use transposed map)
    map_val = tissue_map_t(y, x, z);    
    
    % Convert to colormap index
    map_ind = round(map_val / cmap_sz) + 1;
    map_ind = max(1, min(256, map_ind));
    
    % Calculate vertex index (0-indexed coordinates)
    % vertex_idx = (y_0indexed * nXvoxels + x_0indexed) * 8
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
            
            trisurf(tri, ...
                plot_x, ...
                plot_y, ...
                verts(:, 3), ...
                'LineWidth', 1, ...
                'FaceColor', custom_cmap(map_ind, :), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', alpha);
        catch
            % Skip if delaunay fails (e.g., collinear points)
            continue;
        end
    end
end

end