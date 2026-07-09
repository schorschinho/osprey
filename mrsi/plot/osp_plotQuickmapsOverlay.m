function out = osp_plotQuickmapsOverlay(MRSCont, spec, target, idx_start, idx_end, viridis_map, percentile_clean, target_image, alpha, brainMask, cbar, convention)
%% out = osp_plotQuickmapsOverlay(MRSCont, spec, target, idx_start, idx_end, viridis_map, percentile_clean, target_image, alpha, brainMask, cbar, convention)
%   Creates a figure showing metabolite quickmaps overlaid on MRI images.
%   Handles axial, sagittal, coronal, and oblique acquisitions.
%
%   USAGE:
%       out = osp_plotQuickmapsOverlay(MRSCont, spec, target, idx_start, idx_end, viridis_map, percentile_clean, target_image, alpha, brainMask, cbar, convention)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont         = Osprey data container.
%       spec            = Target spectrum ('raw', etc.)
%       target          = Target metabolite ('tNAA', 'tCr', etc.)
%       idx_start       = index MRSI slice to start (1 at bottom)
%       idx_end         = index MRSI slice to end
%       viridis_map     = use viridis colormap (1) or hot colormap (0)
%       percentile_clean = apply percentile cleanup (0 = off, >0 = percentile value)
%       target_image    = Target image to overlay on:
%                         'MRSIloc', 'T1w_rMRSIloc', 'T1w_rMRSI', 'MRSIloc_rMRSI'
%       alpha           = alpha level of metabolite map overlay (0-1)
%       brainMask       = apply brain mask (0 = off, 1 = apply mask, 2 = show mask only)
%       cbar            = add color bar (0 or 1)
%       convention      = 'radiological' (L on right) or 'neurological' (L on left)
%
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2026-01-21: First version of the code.

%% Parse inputs and set defaults
if nargin < 12
    convention = 'neurological';
    if nargin < 11
        cbar = 1;
        if nargin < 10
            brainMask = 0;
            if nargin < 9
                alpha = 0.66;
                if nargin < 8
                    if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                        target_image = 'T1w_rMRSIloc';
                    else
                        target_image = 'T1w_rMRSI';
                    end
                    if nargin < 7
                        percentile_clean = 0;
                        if nargin < 6
                            viridis_map = 1;
                            if nargin < 5
                                idx_end = MRSCont.raw{1, 1}.nZvoxels;
                                if nargin < 4
                                    idx_start = 1;
                                    if nargin < 3
                                        target = 'tNAA';
                                        if nargin < 2
                                            spec = 'raw';
                                            if nargin < 1
                                                error('ERROR: no input Osprey container specified. Aborting!!');
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Validate convention input
if ~ismember(lower(convention), {'radiological', 'neurological'})
    warning('Invalid convention "%s". Using "radiological" as default.', convention);
    convention = 'radiological';
end
convention = lower(convention);

%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot quickmaps, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot quickmaps, but coregistration has not been performed yet. Run OspreyCoreg first.')
end

%% Validate target exists
if ~isfield(MRSCont.quickMaps, spec)
    available_specs = fieldnames(MRSCont.quickMaps);
    error('Spectrum "%s" not found. Available spectra: %s', spec, strjoin(available_specs, ', '));
end

if ~isfield(MRSCont.quickMaps.(spec), target)
    available_targets = fieldnames(MRSCont.quickMaps.(spec));
    error('Target "%s" not found in spectrum "%s". Available targets: %s', target, spec, strjoin(available_targets, ', '));
end

%% Handle compressed files
[~, ~, T1ext] = fileparts(MRSCont.files_nii{1});
if strcmp(T1ext, '.gz')
    gunzip(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1}, '.gz', '');
end

MRSIlocext = '';
if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    [~, ~, MRSIlocext] = fileparts(MRSCont.files_nii_MRSIloc{1});
    if strcmp(MRSIlocext, '.gz')
        gunzip(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1}, '.gz', '');
    end
end

%% Get gifti vertices
g = gifti(MRSCont.gii_filename_VoxelGrid{1});
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

%% Get quickmap data
quickmap = MRSCont.quickMaps.(spec).(target);

%% Get mask for voxel selection
if MRSCont.flags.didSeg
    voxel_mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
else
    if isfield(MRSCont.opts.MRSI, 'MRSImask')
        voxel_mask = MRSCont.opts.MRSI.MRSImask;
    else
        voxel_mask = ones(MRSCont.raw{1}.nXvoxels, MRSCont.raw{1}.nYvoxels, MRSCont.raw{1}.nZvoxels);
    end
end

%% Apply brain mask to quickmap if requested
switch brainMask
    case 1
        % Apply brain mask
        brain_mask_data = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
        quickmap = quickmap .* brain_mask_data;
    case 2
        % Show brain mask only (for debugging)
        quickmap = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));       
    otherwise
        % No mask applied
        voxel_mask = ones(MRSCont.raw{1}.nXvoxels, MRSCont.raw{1}.nYvoxels, MRSCont.raw{1}.nZvoxels);
end

%% Apply percentile cleanup if requested
if percentile_clean == 1
    quickmap_temp = squeeze(quickmap(:, :, idx_start:idx_end));
    max_val_clean = prctile(quickmap_temp(:), 97);
elseif percentile_clean > 1
    % Use manual max value
    max_val_clean = percentile_clean;
end

%% Calculate colormap range
quickmap_temp = quickmap(:, :, idx_start:idx_end);
if percentile_clean == 1
        valid_values = quickmap_temp(~isnan(quickmap_temp) & plotMap_temp > 0);
        if ~isempty(valid_values)
            max_val_cmap = prctile(valid_values, 97);
        else
            max_val_cmap = 1;
        end
    elseif percentile_clean > 1
        max_val_cmap = percentile_clean;
        else
        max_val_cmap = max(quickmap_temp(:), [], 'all');
end

cmap_sz = max_val_cmap / 255;

%% Setup colormap
if viridis_map
    custom_cmap = viridis(256);
else
    custom_cmap = hot(256);
end

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

%% Display montage of background image
max_val_img = prctile(Coreg_img_display(:), 99.5);

% Get tile dimensions from the display image (after permutation)
tile_height = size(Coreg_img_display, 1);
tile_width = size(Coreg_img_display, 2);

montage(squeeze(Coreg_img_display(:, :, idx_start:idx_end)), ...
    'ThumbnailSize', [tile_height, tile_width], ...
    'Size', [num_rows, slices_per_row], 'DisplayRange', [0 max_val_img]);

% Freeze the grayscale colormap for the background
freezeColors;
hold on;

%% Plot quickmap overlay

osp_plot_quickmap_overlay(quickmap, voxel_mask, vertices_display, ...
    MRSCont.raw{1}.nXvoxels, idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, custom_cmap, cmap_sz, alpha, orientation_info, display_info);

%% Setup colormap and colorbar
colormap(custom_cmap);
clim([0 max_val_cmap]);

if cbar
    cbar_map = colorbar;
    set(cbar_map, 'Color', [1 1 1]);
    cbar_map.Label.String = target;
    cbar_map.Label.Color = [1 1 1];
end

%% Add orientation labels
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

%% Add spectrum/target label
text(10, 15, sprintf('%s - %s', spec, target), ...
    'Rotation', 0, 'Color', 'w', 'FontSize', 15, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');

axis image;
axis off;
hold off;

%% Cleanup compressed files
if strcmp(T1ext, '.gz')
    delete(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1}, '.nii', '.nii.gz');
end

if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    if strcmp(MRSIlocext, '.gz')
        delete(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1}, '.nii', '.nii.gz');
    end
end

end