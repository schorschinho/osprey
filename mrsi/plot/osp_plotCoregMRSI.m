function out = osp_plotCoregMRSI(MRSCont, target_image, addGrid, addOuterMask, addBrainMask, addLipidMask, idx_start, idx_end, convention)
%% out = osp_plotCoregMRSI(MRSCont, target_image, addGrid, addOuterMask, addBrainMask, addLipidMask, idx_start, idx_end, convention)
%   Creates a figure showing the coregistration between MRSI data and MRI
%   images. Handles axial, sagittal, coronal, and oblique acquisitions.
%
%   USAGE:
%       out = osp_plotCoregMRSI(MRSCont, target_image, addGrid, addOuterMask, addBrainMask, addLipidMask, idx_start, idx_end, convention)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont      = Osprey data container.
%       target_image = Target image to overlay on:
%                      'MRSIloc', 'T1w_rMRSIloc', 'T1w_rMRSI', 'MRSIloc_rMRSI'
%       addGrid      = add Grid dots on image (0 or 1)
%       addOuterMask = add outer mask overlay (0 or 1)
%       addBrainMask = add automated brain mask (0 or 1)
%       addLipidMask = add automated lipid mask (0 or 1)
%       idx_start    = index MRSI slice to start (1 at bottom)
%       idx_end      = index MRSI slice to end
%       convention   = 'radiological' (L on right) or 'neurological' (L on left)
%
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

%% Parse inputs and set defaults
if nargin < 9
    convention = 'neurological';
    if nargin < 8
        idx_end = MRSCont.raw{1, 1}.nZvoxels;
        if nargin < 7
            idx_start = 1;
            if nargin < 6
                addLipidMask = 1;
                if nargin < 5
                    addBrainMask = 1;
                    if nargin < 4
                        addOuterMask = 1;
                        if nargin < 3
                            addGrid = 1;
                            if nargin < 2
                                if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                                    target_image = 'T1w_rMRSIloc';
                                else
                                    target_image = 'T1w_rMRSI';
                                end
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

%% Validate convention input
if ~ismember(lower(convention), {'radiological', 'neurological'})
    warning('Invalid convention "%s". Using "radiological" as default.', convention);
    convention = 'radiological';
end
convention = lower(convention);

%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot coregistration, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot coregistration, but coregistration has not been performed yet. Run OspreyCoreg first.')
end

if ~MRSCont.flags.didSeg
    addBrainMask = 0;
    addLipidMask = 0;
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

%% Setup colors
[VoxColors] = cbrewer('qual', 'Set1', 9);

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

%% Display montage
max_val = prctile(Coreg_img_display(:), 99.5);

% Get tile dimensions from the display image (after permutation)
tile_height = size(Coreg_img_display, 1);
tile_width = size(Coreg_img_display, 2);

montage(squeeze(Coreg_img_display(:, :, idx_start:idx_end)), ...
    'ThumbnailSize', [tile_height, tile_width], ...
    'Size', [num_rows, slices_per_row], 'DisplayRange', [0 max_val]);
hold on;

%% Add outer mask overlay
% Note: Do NOT apply osp_apply_mask_orientation_correction here
% The mask dimension transpose is handled inside osp_plot_mask_overlay
if addOuterMask
    osp_plot_mask_overlay(MRSCont.opts.MRSI.outerMask.mask, 0, ...
        vertices_display, MRSCont.raw{1, 1}.nXvoxels, ...
        idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
        VoxColors(1,:), 0.50, orientation_info, display_info);
end

%% Add brain mask overlay
if addBrainMask
    brain_mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
    
    osp_plot_mask_overlay(brain_mask, 1, ...
        vertices_display, MRSCont.raw{1}.nXvoxels, ...
        idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
        VoxColors(3,:), 0.50, orientation_info, display_info);
end

%% Add lipid mask overlay
if addLipidMask
    lip_mask = squeeze(MRSCont.seg.tissue.lip(1,:,:,:));
    
    osp_plot_mask_overlay(lip_mask, 1, ...
        vertices_display, MRSCont.raw{1}.nXvoxels, ...
        idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
        VoxColors(2,:), 0.75, orientation_info, display_info);
end

%% Add grid overlay
if addGrid
    osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, ...
        idx_start, idx_end, slices_per_row, tile_width, tile_height, ...
        orientation_info, display_info);
end

%% Add orientation labels to each slice
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

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