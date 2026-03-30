function osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, idx_start, idx_end, slices_per_row, tile_width, tile_height, orientation_info, display_info)
%% osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, idx_start, idx_end, slices_per_row, tile_width, tile_height, orientation_info, display_info)
%   This function plots MRSI grid vertex markers as an overlay on a
%   montage display of image slices.
%
%   The function iterates through slices in the montage and plots grid
%   vertices that fall within each slice. Vertex positions are transformed
%   from image coordinates to montage coordinates based on the tile layout.
%
%   After display transformations:
%       vertices_display(:,1) corresponds to image rows (plot Y axis)
%       vertices_display(:,2) corresponds to image columns (plot X axis)
%
%   USAGE:
%       osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, idx_start, idx_end, slices_per_row, tile_width, tile_height, orientation_info, display_info);
%
%   INPUTS:
%       vertices_display   = Nx2 matrix of vertex coordinates in display
%                            space after orientation transformations.
%       vertices_MRSI_voxel = Nx3 matrix of vertex coordinates in MRSI
%                            voxel space for slice membership testing.
%       idx_start          = Starting slice index in the montage.
%       idx_end            = Ending slice index in the montage.
%       slices_per_row     = Number of slice tiles per row in the montage.
%       tile_width         = Width of each tile in pixels.
%       tile_height        = Height of each tile in pixels.
%       orientation_info   = Struct from osp_analyze_orientation containing
%                            orientation analysis results.
%       display_info       = Struct from osp_prepare_display_slice containing
%                            display transformation parameters.
%
%   OUTPUTS:
%       None (plots directly to current axes)
%
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    
    grid_color = [254/255, 186/255, 47/255];
    marker_size = 4;
    
    for slice_idx = idx_start:idx_end
        valid_vertices = abs(vertices_MRSI_voxel(:, 3) - slice_idx) < 1;
        
        if ~any(valid_vertices)
            continue;
        end
        
        slice_in_montage = slice_idx - idx_start;
        montage_col = mod(slice_in_montage, slices_per_row);
        montage_row = floor(slice_in_montage / slices_per_row);
        
        montage_x_offset = montage_col * tile_width;
        montage_y_offset = montage_row * tile_height;
        
        % Consistent mapping: dim1 = Y (rows), dim2 = X (cols)
        x_coords = vertices_display(valid_vertices, 2) + montage_x_offset;
        y_coords = vertices_display(valid_vertices, 1) + montage_y_offset;
        
        plot(x_coords, y_coords, '.', 'Color', grid_color, 'MarkerSize', marker_size);
    end
end