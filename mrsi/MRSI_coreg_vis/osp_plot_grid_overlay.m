function osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, ...
    idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
    orientation_info, display_info)
    % Plot grid overlay on montage
    %
    % After all transformations:
    %   vertices_display(:,1) = image rows = plot Y
    %   vertices_display(:,2) = image cols = plot X
    
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