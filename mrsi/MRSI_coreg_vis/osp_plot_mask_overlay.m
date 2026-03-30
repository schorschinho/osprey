function osp_plot_mask_overlay(mask, mask_value, vertices_display, nXvoxels, idx_start, idx_end, slices_per_row, tile_width, tile_height, color, alpha, orientation_info, display_info)
%% osp_plot_mask_overlay
%   Plots mask overlay on montage.
%
%   The masks from MRSCont.seg.tissue are in MRSI voxel space [nX, nY, nZ].
%   The vertex indexing uses: vertex_idx = (y * nXvoxels + x) * 8
%   So we need to swap dimensions 1 and 2 to match.
%
%   USAGE:
%       osp_plot_mask_overlay(mask, mask_value, vertices_display, nXvoxels, ...
%           idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
%           color, alpha, orientation_info, display_info)
%
%   INPUTS:
%       mask            = 3D mask in MRSI voxel space [nX, nY, nZ]
%       mask_value      = 0 to plot where mask==0, 1 to plot where mask==1
%       vertices_display = transformed vertex coordinates for display
%       nXvoxels        = number of X voxels in MRSI grid
%       idx_start       = start slice index
%       idx_end         = end slice index
%       slices_per_row  = number of slices per row in montage
%       tile_width      = width of each tile in pixels
%       tile_height     = height of each tile in pixels
%       color           = RGB color for overlay
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
%% Transpose mask to match vertex ordering
% Masks are stored as [nX, nY, nZ] but vertex indexing assumes [nY, nX, nZ]
% vertex_idx = (y * nXvoxels + x) * 8, where y is outer loop, x is inner
mask_t = permute(mask, [2, 1, 3]);  % Now [nY, nX, nZ]

%% Find voxels to plot
if mask_value == 0
    linear_indices = find(mask_t == 0);
else
    linear_indices = find(mask_t);
end

if isempty(linear_indices)
    return;
end

%% Get voxel indices (now in [nY, nX, nZ] order after permute)
[yIdx, xIdx, zIdx] = ind2sub(size(mask_t), linear_indices);

%% Filter to slice range
valid_slices = (zIdx >= idx_start) & (zIdx <= idx_end);
xIdx = xIdx(valid_slices);
yIdx = yIdx(valid_slices);
zIdx = zIdx(valid_slices);

if isempty(xIdx)
    return;
end

%% Plot each voxel
for vox = 1:length(xIdx)
    % Current voxel indices (1-indexed)
    x = xIdx(vox);  % X index (1 to nXvoxels)
    y = yIdx(vox);  % Y index (1 to nYvoxels)
    z = zIdx(vox);  % Z index (slice)
    
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
                'FaceColor', color, ...
                'EdgeColor', 'none', ...
                'FaceAlpha', alpha);
        catch
            % Skip if delaunay fails (e.g., collinear points)
            continue;
        end
    end
end

end