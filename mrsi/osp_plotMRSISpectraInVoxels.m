function out = osp_plotMRSISpectraInVoxels(MRSCont, target_image, kk, which, idx_start, idx_end, convention, options)
%% out = osp_plotMRSISpectraInVoxels(MRSCont, target_image, kk, which, idx_start, idx_end, convention, options)
%   Creates a figure showing MRSI spectra plotted inside each voxel region
%   overlaid on coregistered MRI images.
%
%   USAGE:
%       out = osp_plotMRSISpectraInVoxels(MRSCont, target_image, kk, which, idx_start, idx_end, convention, options)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   ARGUMENTS:
%       MRSCont      = Osprey data container.
%       target_image = Target image to overlay on:
%                      'MRSIloc', 'T1w_rMRSIloc', 'T1w_rMRSI', 'MRSIloc_rMRSI'
%       kk           = Index of the dataset (default: 1)
%       which        = Which data to plot: 'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum' (default: 'A')
%       idx_start    = index MRSI slice to start (1 at bottom)
%       idx_end      = index MRSI slice to end
%       convention   = 'radiological' (L on right) or 'neurological' (L on left)
%       options      = struct with additional options:
%                      .ppm_range     = [min max] ppm range to display (default: [0.5 4.2])
%                      .spec_color    = spectrum line color (default: 'g')
%                      .line_width    = spectrum line width (default: 0.5)
%                      .show_grid     = show voxel grid outlines (default: true)
%                      .grid_color    = grid outline color (default: [0.5 0.5 0.5])
%                      .fill_factor   = how much of voxel to fill with spectrum (default: 0.8)
%                      .threshold     = SNR or amplitude threshold to show spectrum (default: 0)
%                      .show_background = show MRI background (default: true)
%                      .invert_spec   = invert spectrum display (default: false)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025)
%       hzoelln2@jhmi.edu

%% Parse inputs and set defaults
if nargin < 8 || isempty(options)
    options = struct();
end

if nargin < 7 || isempty(convention)
    convention = 'neurological';
end

if nargin < 6 || isempty(idx_end)
    idx_end = MRSCont.raw{1, 1}.nZvoxels;
end

if nargin < 5 || isempty(idx_start)
    idx_start = 1;
end

if nargin < 4 || isempty(which)
    which = 'A';
end

if nargin < 3 || isempty(kk)
    kk = 1;
end

if nargin < 2 || isempty(target_image)
    if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
        target_image = 'T1w_rMRSIloc';
    else
        target_image = 'T1w_rMRSI';
    end
end

if nargin < 1
    error('ERROR: no input Osprey container specified. Aborting!!');
end

%% Set default options
default_options = struct(...
    'ppm_range', [0.5 4.2], ...
    'spec_color', [254/255, 186/255, 47/255], ...
    'line_width', 0.75, ...
    'show_grid', true, ...
    'grid_color', [0.5 0.5 0.5], ...
    'fill_factor', 0.8, ...
    'threshold', 0, ...
    'show_background', true, ...
    'invert_spec', true, ...
    'use_real', true);

% Merge user options with defaults
opt_fields = fieldnames(default_options);
for i = 1:length(opt_fields)
    if ~isfield(options, opt_fields{i})
        options.(opt_fields{i}) = default_options.(opt_fields{i});
    end
end

%% Validate convention input
if ~ismember(lower(convention), {'radiological', 'neurological'})
    warning('Invalid convention "%s". Using "neurological" as default.', convention);
    convention = 'neurological';
end
convention = lower(convention);

%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot spectra, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot coregistration, but coregistration has not been performed yet. Run OspreyCoreg first.')
end

%% Handle compressed files
[~, ~, T1ext] = fileparts(MRSCont.files_nii{kk});
if strcmp(T1ext, '.gz')
    gunzip(MRSCont.files_nii{kk});
    MRSCont.files_nii{kk} = strrep(MRSCont.files_nii{kk}, '.gz', '');
end

MRSIlocext = '';
if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{kk})
    [~, ~, MRSIlocext] = fileparts(MRSCont.files_nii_MRSIloc{kk});
    if strcmp(MRSIlocext, '.gz')
        gunzip(MRSCont.files_nii_MRSIloc{kk});
        MRSCont.files_nii_MRSIloc{kk} = strrep(MRSCont.files_nii_MRSIloc{kk}, '.gz', '');
    end
end

%% Get MRSI data dimensions
nXvoxels = MRSCont.raw{kk}.nXvoxels;
nYvoxels = MRSCont.raw{kk}.nYvoxels;
nZvoxels = MRSCont.raw{kk}.nZvoxels;

%% Get spectral data
[spectra, ppm] = osp_get_mrsi_spectra(MRSCont, kk, which);

%% Get gifti vertices
g = gifti(MRSCont.gii_filename_VoxelGrid{kk});
vertices = g.vertices;

%% Load and prepare images based on target_image type
switch target_image
    case 'MRSIloc'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices);
        
    case 'T1w_rMRSIloc'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices);
        
    case 'T1w_rMRSI'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSI(MRSCont, vertices);
        
    case 'MRSIloc_rMRSI'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc_rMRSI(MRSCont, vertices);
        
    otherwise
        error('Unknown target_image type: %s', target_image);
end

%% Analyze orientation
[orientation_info] = osp_analyze_orientation(AffineMat);

%% Prepare image for display (with permutation and flips)
[Coreg_img_display, display_info] = osp_prepare_image_for_display(Coreg_img, orientation_info, convention);

%% Transform vertices to match displayed image
[vertices_display] = osp_transform_vertices_for_display(vertices_voxel, Coreg_img, orientation_info, display_info, MRSI_vol);

% Also transform for MRSI voxel space (needed for slice selection)
vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
vertices_MRSI_voxel = (MRSI_vol.mat \ vertices_homogeneous')';
vertices_MRSI_voxel = vertices_MRSI_voxel(:, 1:3);


%% Setup figure
out = figure;
set(out, 'Color', [0 0 0]);

%% Calculate montage layout
num_slices = idx_end - idx_start + 1;
if num_slices < 5
    slices_per_row = num_slices;
else
    slices_per_row = 5;
end
num_rows = ceil(num_slices / slices_per_row);

%% Display montage background
max_val = prctile(Coreg_img_display(:), 99.5);

% Get tile dimensions from the display image (after permutation)
tile_height = size(Coreg_img_display, 1);
tile_width = size(Coreg_img_display, 2);

if options.show_background
    montage(squeeze(Coreg_img_display(:, :, idx_start:idx_end)), ...
        'ThumbnailSize', [tile_height, tile_width], ...
        'Size', [num_rows, slices_per_row], 'DisplayRange', [0 max_val]);
else
    % Create black background
    montage(zeros(size(squeeze(Coreg_img_display(:, :, idx_start:idx_end)))), ...
        'ThumbnailSize', [tile_height, tile_width], ...
        'Size', [num_rows, slices_per_row]);
end
hold on;
%%


%% Plot spectra in each voxel
osp_plot_spectra_in_voxels(spectra, ppm, vertices_display, vertices_MRSI_voxel, ...
    nXvoxels, nYvoxels, nZvoxels, ...
    idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
    orientation_info, display_info, options);

%% Add grid overlay if requested
if options.show_grid
    osp_plot_voxel_grid_overlay(vertices_display, vertices_MRSI_voxel, ...
        nXvoxels, nYvoxels, nZvoxels, ...
        idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
        orientation_info, display_info, options.grid_color);
end

%% Add orientation labels to each slice
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

axis image;
axis off;
hold off;

%% Cleanup compressed files
if strcmp(T1ext, '.gz')
    delete(MRSCont.files_nii{kk});
    MRSCont.files_nii{kk} = strrep(MRSCont.files_nii{kk}, '.nii', '.nii.gz');
end

if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{kk})
    if strcmp(MRSIlocext, '.gz')
        delete(MRSCont.files_nii_MRSIloc{kk});
        MRSCont.files_nii_MRSIloc{kk} = strrep(MRSCont.files_nii_MRSIloc{kk}, '.nii', '.nii.gz');
    end
end

end

%% ========================================================================
%% Helper Functions
%% ========================================================================

function [spectra, ppm] = osp_get_mrsi_spectra(MRSCont, kk, which)
%% osp_get_mrsi_spectra - Extract spectral data from MRSCont
%
%   Returns spectra as 4D array: [nPoints, nXvoxels, nYvoxels, nZvoxels]

% Get the appropriate data based on 'which' parameter
switch upper(which)
    case 'A'
        if isfield(MRSCont, 'processed') && isfield(MRSCont.processed, 'A')
            data = MRSCont.processed.A{kk};
        else
            data = MRSCont.raw{kk};
        end
    case 'AFID'
        data = MRSCont.processed.AFID{kk};
    case 'B'
        data = MRSCont.processed.B{kk};
    case 'C'
        data = MRSCont.processed.C{kk};
    case 'D'
        data = MRSCont.processed.D{kk};
    case 'DIFF1'
        data = MRSCont.processed.diff1{kk};
    case 'DIFF2'
        data = MRSCont.processed.diff2{kk};
    case 'SUM'
        data = MRSCont.processed.sum{kk};
    otherwise
        error('Unknown data type: %s', which);
end

% Extract ppm axis
ppm = data.ppm;

% Extract spectra - handle different data structures
if isfield(data, 'specs')
    spectra = data.specs;
elseif isfield(data, 'fids')
    % Convert FIDs to spectra
    spectra = fftshift(fft(data.fids, [], 1), 1);
else
    error('Cannot find spectral data in MRSCont structure');
end

% Ensure correct dimensions [nPoints, nX, nY, nZ]
dims = size(spectra);
if length(dims) == 2
    % Reshape from [nPoints, nVoxels] to [nPoints, nX, nY, nZ]
    nPoints = dims(1);
    nVoxels = dims(2);
    nX = data.nXvoxels;
    nY = data.nYvoxels;
    nZ = data.nZvoxels;
    spectra = reshape(spectra, [nPoints, nX, nY, nZ]);
end

end

function osp_plot_spectra_in_voxels(spectra, ppm, vertices_display, vertices_MRSI_voxel, ...
    nXvoxels, nYvoxels, nZvoxels, ...
    idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
    orientation_info, display_info, options)
%% osp_plot_spectra_in_voxels - Plot spectra inside voxel boundaries
%
%   MRSI coordinate to display mapping:
%   - MRSI X (1->34): top to bottom in display
%   - MRSI Y (1->44): right to left in display
%   - MRSI Z: slice index
%
%   Spectra indexing: spectra(:, MRSI_X, MRSI_Y, MRSI_Z)

% Get ppm indices for display range
ppm_idx = ppm >= options.ppm_range(1) & ppm <= options.ppm_range(2);
nPpmPoints = sum(ppm_idx);

% Number of vertices per voxel (8 for hexahedron)
vertices_per_voxel = 8;
total_voxels = size(vertices_display, 1) / vertices_per_voxel;

% Loop through each voxel directly using vertices
for voxel_idx = 1:total_voxels
    
    % Get vertex indices for this voxel
    vert_start = (voxel_idx - 1) * vertices_per_voxel + 1;
    vert_end = voxel_idx * vertices_per_voxel;
    
    % Get display vertices for this voxel
    voxel_vertices = vertices_display(vert_start:vert_end, :);
    
    % Get MRSI voxel space vertices for this voxel
    voxel_mrsi_coords = vertices_MRSI_voxel(vert_start:vert_end, :);
    
    % Skip if vertices contain NaN
    if any(isnan(voxel_vertices(:))) || any(isnan(voxel_mrsi_coords(:)))
        continue;
    end
    
    % Calculate center in MRSI voxel space
    mrsi_center = mean(voxel_mrsi_coords, 1);
    
    % Get spectrum indices directly from MRSI coordinates
    mrsi_x = round(mrsi_center(1));
    mrsi_y = round(mrsi_center(2));
    mrsi_z = round(mrsi_center(3));
    
    % Check if this slice is in our display range
    if mrsi_z < idx_start || mrsi_z > idx_end
        continue;
    end
    
    % Validate indices
    if mrsi_x < 1 || mrsi_x > size(spectra, 2) || ...
       mrsi_y < 1 || mrsi_y > size(spectra, 3) || ...
       mrsi_z < 1 || mrsi_z > size(spectra, 4)
        continue;
    end
    
    % Calculate montage position
    montage_idx = mrsi_z - idx_start;
    montage_row = floor(montage_idx / slices_per_row);
    montage_col = mod(montage_idx, slices_per_row);
    x_offset = montage_col * tile_width;
    y_offset = montage_row * tile_height;
    
    % Get spectrum for this voxel
    % spectra is indexed as: spectra(:, MRSI_X, MRSI_Y, MRSI_Z)
    spectrum = squeeze(spectra(:, mrsi_x, mrsi_y, mrsi_z));
    
    % Use real part if specified
    if options.use_real
        spectrum = real(spectrum);
    else
        spectrum = abs(spectrum);
    end
    
    % Extract ppm range
    spectrum = spectrum(ppm_idx);
    
    % Skip if spectrum is all zeros or below threshold
    if all(spectrum == 0) || max(abs(spectrum)) < options.threshold
        continue;
    end
    
    % Normalize spectrum if requested
    % if options.normalize
        spectrum = spectrum - min(spectrum);
        max_spec = max(spectrum);
        if max_spec > 0
            spectrum = spectrum / max_spec;
        end
    % end
    
    % Invert if requested
    % if options.invert_spec
        spectrum = 1 - spectrum;
    % end
    
    % Get voxel bounding box in display coordinates
    vox_x = voxel_vertices(:, 1);
    vox_y = voxel_vertices(:, 2);
    
    % Get 2D convex hull of the voxel projection
    try
        k = convhull(vox_x, vox_y);
        hull_x = vox_x(k);
        hull_y = vox_y(k);
    catch
        hull_x = [min(vox_x); max(vox_x); max(vox_x); min(vox_x); min(vox_x)];
        hull_y = [min(vox_y); min(vox_y); max(vox_y); max(vox_y); min(vox_y)];
    end
    
    % Calculate bounding box
    vox_x_min = min(hull_x);
    vox_x_max = max(hull_x);
    vox_y_min = min(hull_y);
    vox_y_max = max(hull_y);
    
    vox_width = vox_x_max - vox_x_min;
    vox_height = vox_y_max - vox_y_min;
    
    % Apply fill factor
    center_x = (vox_x_min + vox_x_max) / 2;
    center_y = (vox_y_min + vox_y_max) / 2;
    
    plot_width = vox_width * options.fill_factor;
    plot_height = vox_height * options.fill_factor;
    
    plot_x_min = center_x - plot_width / 2;
    plot_x_max = center_x + plot_width / 2;
    plot_y_min = center_y - plot_height / 2;
    
    % Map spectrum to voxel coordinates
    % Plot spectrum horizontally within the voxel
    spec_plot_x = linspace(plot_x_max, plot_x_min, nPpmPoints);
    spec_plot_y = plot_y_min + spectrum * plot_height;
    
    % Add tile offset
    spec_plot_x = spec_plot_x + x_offset;
    spec_plot_y = spec_plot_y + y_offset;
    
    % Plot spectrum
    plot(spec_plot_x, spec_plot_y, 'Color', options.spec_color, ...
        'LineWidth', options.line_width);
    
end

end

function osp_plot_voxel_grid_overlay(vertices_display, vertices_MRSI_voxel, ...
    nXvoxels, nYvoxels, nZvoxels, ...
    idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
    orientation_info, display_info, grid_color)
%% osp_plot_voxel_grid_overlay - Plot voxel grid outlines
%
%   This function draws the outline of each voxel

% Total number of voxels per slice
voxels_per_slice = nXvoxels * nYvoxels;

% Number of vertices per voxel
vertices_per_voxel = 8;

% Loop through each MRSI slice
for slice_idx = idx_start:idx_end
    
    % Calculate position in montage
    montage_idx = slice_idx - idx_start;
    montage_row = floor(montage_idx / slices_per_row);
    montage_col = mod(montage_idx, slices_per_row);
    
    % Calculate offset for this tile in the montage
    x_offset = montage_col * tile_width;
    y_offset = montage_row * tile_height;
    
    % Loop through each voxel in this slice
    for y_vox = 1:nYvoxels
        for x_vox = 1:nXvoxels
            
            % Calculate linear voxel index
            voxel_idx = (slice_idx - 1) * voxels_per_slice + ...
                        (y_vox - 1) * nXvoxels + x_vox;
            
            % Get vertex indices for this voxel
            vert_start = (voxel_idx - 1) * vertices_per_voxel + 1;
            vert_end = voxel_idx * vertices_per_voxel;
            
            if vert_end > size(vertices_display, 1)
                continue;
            end
            
            % Get vertices for this voxel
            voxel_vertices = vertices_display(vert_start:vert_end, :);
            
            % Skip if vertices contain NaN
            if any(isnan(voxel_vertices(:)))
                continue;
            end
            
            % Get 2D projection
            vox_x = voxel_vertices(:, 1);
            vox_y = voxel_vertices(:, 2);
            
            % Get 2D convex hull
            try
                k = convhull(vox_x, vox_y);
                hull_x = vox_x(k) + x_offset;
                hull_y = vox_y(k) + y_offset;
                
                plot(hull_x, hull_y, 'Color', grid_color, 'LineWidth', 0.5);
            catch
                % If convhull fails, draw bounding box
                bx = [min(vox_x); max(vox_x); max(vox_x); min(vox_x); min(vox_x)] + x_offset;
                by = [min(vox_y); min(vox_y); max(vox_y); max(vox_y); min(vox_y)] + y_offset;
                plot(bx, by, 'Color', grid_color, 'LineWidth', 0.5);
            end
            
        end
    end
end

end