function [out_x_cord, out_y_cord] = translateCoordinatesBetweenTools(x_cord, y_cord, nXvoxels, nYvoxels, direction)
    % Translates between FSL and Osprey coordinate systems with 90° rotation
    %
    % Inputs:
    %   x_cord, y_cord: input coordinates
    %   nXvoxels, nYvoxels: ORIGINAL matrix dimensions (always in Osprey orientation)
    %   direction: 'FSLtoOsprey' or 'OspreytoFSL'
    %
    % Outputs:
    %   out_x_cord, out_y_cord: translated coordinates
    %
    % Note: FSL uses 0-based indexing with origin at bottom-right
    %       Osprey uses 1-based indexing with origin at top-left, rotated 90° CW
    
    switch direction
        case 'FSLtoOsprey'
            % FSL dimensions are swapped due to rotation: (nYvoxels, nXvoxels)
            % Inverse of the OspreytoFSL transformation
            out_x_cord = nXvoxels - x_cord;
            out_y_cord = y_cord + 1;
            
        case 'OspreytoFSL'
            % Osprey to FSL with 90° CW rotation
            out_y_cord = y_cord - 1;  % -1 for 1-based to 0-based conversion
            out_x_cord = nXvoxels - x_cord;
            
        otherwise
            error('Direction must be ''FSLtoOsprey'' or ''OspreytoFSL''');
    end
end