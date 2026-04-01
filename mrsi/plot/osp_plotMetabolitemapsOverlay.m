function out = osp_plotMetabolitemapsOverlay(MRSCont, quantification, metabolite, idx_start, idx_end, viridis_map, percentile_clean, clip, target_image, alpha, cbar, convention)
%% out = osp_plotMetabolitemapsOverlay(MRSCont, quantification, metabolite, idx_start, idx_end, viridis_map, percentile_clean, clip, target_image, alpha, cbar, convention)
%   Creates a figure showing metabolite maps overlaid on MRI images.
%   Handles axial, sagittal, coronal, and oblique acquisitions.
%
%   USAGE:
%       out = osp_plotMetabolitemapsOverlay(MRSCont, quantification, metabolite, idx_start, idx_end, viridis_map, percentile_clean, clip, target_image, alpha, cbar, convention)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont         = Osprey data container.
%       quantification  = Which quantification to plot ('amplitudes', 'water', 'GlobalQC', etc.)
%       metabolite      = Target metabolite ('tNAA_Acetyl_only', 'tCr', etc.)
%       idx_start       = index MRSI slice to start (1 at bottom)
%       idx_end         = index MRSI slice to end
%       viridis_map     = use viridis colormap (1) or hot colormap (0)
%       percentile_clean = apply percentile cleanup (0 = off, 1 = 97th, percentile, >1 = manual max value)\
%       clip            = clip colormap (1 = clip, 0 = remove)
%       target_image    = Target image to overlay on:
%                         'MRSIloc', 'T1w_rMRSIloc', 'T1w_rMRSI', 'MRSIloc_rMRSI'
%       alpha           = alpha level of metabolite map overlay (0-1)
%       cbar            = add color bar (0 or 1)
%       convention      = 'radiological' (L on right) or 'neurological' (L on left)
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
            alpha = 0.66;
            if nargin < 9
                if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                    target_image = 'T1w_rMRSIloc';
                else
                    target_image = 'T1w_rMRSI';
                end
                if nargin < 8
                    clip = 1;
                    if nargin < 7
                        percentile_clean = 0;
                        if nargin < 6
                            viridis_map = 1;
                            if nargin < 5
                                idx_end = MRSCont.raw{1, 1}.nZvoxels;
                                if nargin < 4
                                    idx_start = 1;
                                    if nargin < 3
                                        metabolite = 'tNAA_Acetyl_only';
                                        if nargin < 2
                                            quantification = 'amplitudes';
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
    warning('Invalid convention "%s". Using "neurological" as default.', convention);
    convention = 'neurological';
end
convention = lower(convention);

%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot metabolite maps, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot metabolite maps, but coregistration has not been performed yet. Run OspreyCoreg first.')
end

%% Validate quantification and metabolite exist
if ~isfield(MRSCont.quantify, quantification)
    available_quants = fieldnames(MRSCont.quantify);
    error('Quantification "%s" not found. Available quantifications: %s', quantification, strjoin(available_quants, ', '));
end

if ~strcmp(quantification, 'water') && ~strcmp(quantification, 'GlobalQC')
    if ~isfield(MRSCont.quantify.(quantification), metabolite)
        available_metabolites = fieldnames(MRSCont.quantify.(quantification));
        error('Metabolite "%s" not found in quantification "%s". Available metabolites: %s', metabolite, quantification, strjoin(available_metabolites, ', '));
    end
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

%% Get metabolite map data
if ~strcmp(quantification, 'water') && ~strcmp(quantification, 'GlobalQC')
    plotMap = MRSCont.quantify.(quantification).(metabolite);
else
    plotMap = MRSCont.quantify.(quantification);
end

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

%% Determine if this is a QC plot
isQC = contains(quantification, 'QC') && ~contains(quantification, 'QCfilt');

%% Apply percentile cleanup if requested (only for non-QC plots)
if ~isQC
    if percentile_clean == 1
        plotMap_temp = squeeze(plotMap(:, :, idx_start:idx_end));
        max_val_clean = prctile(plotMap_temp(:), 97);
        if ~clip
            plotMap(plotMap > max_val_clean) = NaN;
        end
    elseif percentile_clean > 1
        % Use manual max value
        max_val_clean = percentile_clean;
        if ~clip
            plotMap(plotMap > max_val_clean) = NaN;
        end
    end
end

%% Calculate colormap range
plotMap_temp = plotMap(:, :, idx_start:idx_end);
if ~isQC
    if percentile_clean == 1
        valid_values = plotMap_temp(~isnan(plotMap_temp) & plotMap_temp > 0);
        if ~isempty(valid_values)
            max_val_cmap = prctile(valid_values, 97);
        else
            max_val_cmap = 1;
        end
    elseif percentile_clean > 1
        max_val_cmap = percentile_clean;
        else
        max_val = max(plotMapTemp(:), [], 'all');
    end
else
    max_val_cmap = max(plotMap_temp(:), [], 'all');
end
cmap_sz = max_val_cmap / 255;

%% Setup colormap
if isQC
    % Special colormaps for QC plots
    if strcmp(quantification, 'GlobalQC')
        custom_cmap = [0 0 0;
                       1 0 0;
                       0 0 1;
                       0 1 0];
    else
        custom_cmap = [0 0 0;
                       1 0 0;
                       0 0 1;
                       1 0 1;
                       0 1 1;
                       0 1 0];
    end
else
    % Standard colormaps for metabolite maps
    if viridis_map
        custom_cmap = viridis(256);
    else
        custom_cmap = hot(256);
    end
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

%% Plot metabolite map overlay
osp_plot_quickmap_overlay(plotMap, voxel_mask, vertices_display, ...
    MRSCont.raw{1}.nXvoxels, idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, custom_cmap, cmap_sz, alpha, orientation_info, display_info, isQC);

%% Setup colormap and colorbar
colormap(custom_cmap);
clim([0 max_val_cmap]);

if cbar
    cbar_map = colorbar;
    set(cbar_map, 'Color', [1 1 1]);

    if isQC
        % QC-specific colorbar labels
        cbar_map.Label.String = 'Applied QC filter';
        cbar_map.Label.Color = [1 1 1];
        if strcmp(quantification, 'GlobalQC')
            cbar_map.Ticks = [0 1 2 3];
            cbar_map.TickLabels = {'brain mask', 'FWHM', 'SNR', 'QC passed'};
        else
            cbar_map.Ticks = [0 1 2 3 4 5];
            cbar_map.TickLabels = {'brain mask', 'FWHM', 'SNR', 'CRLB', 'percentile', 'QC passed'};
        end
    else
        % Standard metabolite colorbar labels
        if strcmp(quantification, 'water')
            cbar_map.Label.String = 'amplitude water';
        else
            cbar_map.Label.String = quantification;
        end
        cbar_map.Label.Color = [1 1 1];
    end
end

%% Add orientation labels
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

%% Add quantification/metabolite label
if ~isQC && ~strcmp(quantification, 'water')
    text(10, 15, sprintf('%s - %s', quantification, metabolite), ...
        'Rotation', 0, 'Color', 'w', 'FontSize', 15, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

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
