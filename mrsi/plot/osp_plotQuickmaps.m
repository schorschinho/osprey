%% out = osp_plotQuickmaps(MRSCont, spec, target, idx_start, idx_end, viridis_map, interpolated, brainMask, convention, cbar)
%   This function creates a figure showing metabolite quickmaps from an
%   MRSI scan as a montage display without anatomical overlay.
%
%   The function handles axial, sagittal, coronal, and oblique acquisitions
%   by automatically detecting the image orientation and applying appropriate
%   transformations for correct anatomical display.
%
%
%   INPUTS:
%       MRSCont      = Osprey MRS data container.
%       spec         = Target spectrum string ('raw', 'A', etc.)
%       target       = Target metabolite string ('tNAA', 'tCr', 'Glx', etc.)
%       idx_start    = Starting MRSI slice index (optional, default: 1)
%       idx_end      = Ending MRSI slice index (optional, default: nZvoxels)
%       viridis_map  = Use viridis colormap (1) or gray colormap (0)
%                      (optional, default: 1)
%       interpolated = Use interpolated maps (0 = off, >0 = use interpolated)
%                      (optional, default: 0)
%       brainMask    = Brain mask application (optional, default: 1):
%                      0 = no mask
%                      1 = apply brain mask
%                      2 = show mask only (for debugging)
%       convention   = Display convention (optional, default: 'neurological'):
%                      'radiological' - Left on right side of image
%                      'neurological' - Left on left side of image
%       cbar         = Display colorbar (0 or 1) (optional, default: 1)
%
%   OUTPUTS:
%       out          = MATLAB figure handle.
%
%   PREREQUISITES:
%       - OspreyLoad must have been run for 'raw' or 'raw_w' spectra
%       - OspreyProcess must have been run for processed spectra
%       - OspreySeg must have been run for brainMask option
%       - OspreyCoreg must have been run for orientation detection
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.
%%

function out = osp_plotQuickmaps(MRSCont, spec, target, idx_start, idx_end, viridis_map, interpolated, brainMask, convention, cbar)

%% Validate prerequisites
if ~MRSCont.flags.didLoadData && (strcmp(spec,'raw') || strcmp(spec,'raw_w'))
    error('Trying to plot quick maps, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didProcess && (strcmp(spec,'A'))
    error('Trying to plot quick maps, but data has not been processed yet. Run OspreyProcess first.')
end

%% Fall back to defaults if not provided
if nargin < 10
    cbar = 1;
    if nargin < 9
        convention = 'neurological';
        if nargin < 8
            brainMask = 1;    
            if nargin < 7
                interpolated = 0;
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

%% Validate convention input
if ~ismember(lower(convention), {'radiological', 'neurological'})
    warning('Invalid convention "%s". Using "neurological" as default.', convention);
    convention = 'neurological';
end
convention = lower(convention);

%% Handle interpolated and brainMask conflicts
if interpolated
    brainMask = 0;
end

if brainMask && ~MRSCont.flags.didSeg
    warning('Trying to apply a brain mask, but data has not been segmented yet. Run OspreySeg first.')
    brainMask = 0;
end

%% Get orientation information
% Try to get affine matrix from available sources
AffineMat = [];
if isfield(MRSCont, 'coreg') && isfield(MRSCont.coreg, 'vol_image') && ...
        ~isempty(MRSCont.coreg.vol_image) && isfield(MRSCont.coreg.vol_image{1}, 'mat')
    AffineMat = MRSCont.coreg.vol_image{1}.mat;
elseif isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc) && ...
        ~isempty(MRSCont.files_nii_MRSIloc{1})
    % Handle compressed files
    [~, ~, ext] = fileparts(MRSCont.files_nii_MRSIloc{1});
    temp_file = MRSCont.files_nii_MRSIloc{1};
    if strcmp(ext, '.gz')
        gunzip(temp_file);
        temp_file = strrep(temp_file, '.gz', '');
        temp_vol = spm_vol(temp_file);
        AffineMat = temp_vol.mat;
        delete(temp_file);
    else
        temp_vol = spm_vol(temp_file);
        AffineMat = temp_vol.mat;
    end
elseif isfield(MRSCont, 'files_nii') && ~isempty(MRSCont.files_nii) && ...
        ~isempty(MRSCont.files_nii{1})
    % Handle compressed files
    [~, ~, ext] = fileparts(MRSCont.files_nii{1});
    temp_file = MRSCont.files_nii{1};
    if strcmp(ext, '.gz')
        gunzip(temp_file);
        temp_file = strrep(temp_file, '.gz', '');
        temp_vol = spm_vol(temp_file);
        AffineMat = temp_vol.mat;
        delete(temp_file);
    else
        temp_vol = spm_vol(temp_file);
        AffineMat = temp_vol.mat;
    end
end

%% Analyze orientation if affine matrix available
has_orientation_info = ~isempty(AffineMat);
if has_orientation_info
    orientation_info = osp_analyze_orientation(AffineMat);
else
    % Create default orientation info for axial
    orientation_info.slice_orientation = 'axial';
    orientation_info.standard_horiz_world = 1;  % L-R
    orientation_info.standard_vert_world = 2;   % A-P
    orientation_info.img_dim_to_world = [1, 2, 3];
    orientation_info.world_to_img_dim = [1, 2, 3];
    orientation_info.dir_cos = eye(3);
    warning('No affine matrix available. Assuming axial orientation.');
end

%% Setup figure
if ~MRSCont.flags.isGUI
    out = figure;   
else
    out = figure('Visible','off');
end

%% Get colormap
if ~viridis_map
    set(out, 'Color', [0 0 0]); 
else
    vir = viridis;
    set(out, 'Color', vir(1,:)); 
end

%% Get output map
if interpolated
    plotMap = MRSCont.quickMapsInt.(spec).(target);
else
    plotMap = MRSCont.quickMaps.(spec).(target);
    switch brainMask
        case 1
            mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
            plotMap = plotMap .* mask;
        case 2
            mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
            plotMap = mask; % For debugging it might be useful to plot the brain mask only
        otherwise
            % No mask applied
    end
end

%% Prepare map for display with correct orientation
% Create a pseudo display_info structure for the quickmap
display_info.convention = convention;
display_info.flipped_dims = [];
display_info.permute_order = [1, 2, 3];
display_info.permuted = false;

% Get the standard arrangement for this orientation
std_horiz = orientation_info.standard_horiz_world;
std_vert = orientation_info.standard_vert_world;
world_to_img_dim = orientation_info.world_to_img_dim;

% Find which image dimensions currently correspond to these world axes
current_horiz_dim = world_to_img_dim(std_horiz);
current_vert_dim = world_to_img_dim(std_vert);
slice_dim = 3;

% Build permutation order
permute_order = [current_vert_dim, current_horiz_dim, slice_dim];

% Check if permutation is needed and valid
needs_permutation = ~isequal(permute_order, [1, 2, 3]);
is_valid_permutation = length(unique(permute_order)) == 3;

if needs_permutation && is_valid_permutation
    display_info.permuted = true;
    display_info.permute_order = permute_order;
    plotMap_display = permute(plotMap, permute_order);
else
    plotMap_display = plotMap;
end

% Store the world axes for the display dimensions
display_info.display_vert_world = std_vert;
display_info.display_horiz_world = std_horiz;

% Get direction signs
if display_info.permuted
    perm = display_info.permute_order;
    vert_img_dim_original = perm(1);
    horiz_img_dim_original = perm(2);
else
    vert_img_dim_original = 1;
    horiz_img_dim_original = 2;
end

dir_cos = orientation_info.dir_cos;
vert_sign = sign(dir_cos(std_vert, vert_img_dim_original));
horiz_sign = sign(dir_cos(std_horiz, horiz_img_dim_original));

display_info.vert_sign_original = vert_sign;
display_info.horiz_sign_original = horiz_sign;

%% Apply flips based on orientation and convention
switch orientation_info.slice_orientation
    case 'axial'
        % Vertical = A-P: A should be at top
        if vert_sign > 0
            plotMap_display = flip(plotMap_display, 1);
            display_info.flipped_dims = [display_info.flipped_dims, 1];
        end
        
        % Horizontal = L-R: depends on convention
        if strcmp(convention, 'radiological')
            if horiz_sign > 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        else  % neurological
            if horiz_sign < 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        end
        
    case 'sagittal'
        % Vertical = S-I: S should be at top
        if vert_sign > 0
            plotMap_display = flip(plotMap_display, 1);
            display_info.flipped_dims = [display_info.flipped_dims, 1];
        end
        
        % Horizontal = A-P: A should be at left
        if horiz_sign > 0
            plotMap_display = flip(plotMap_display, 2);
            display_info.flipped_dims = [display_info.flipped_dims, 2];
        end
        
    case 'coronal'
        % Vertical = S-I: S should be at top
        if vert_sign > 0
            plotMap_display = flip(plotMap_display, 1);
            display_info.flipped_dims = [display_info.flipped_dims, 1];
        end
        
        % Horizontal = L-R: depends on convention
        if strcmp(convention, 'radiological')
            if horiz_sign > 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        else
            if horiz_sign < 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        end
        
    case 'oblique'
        % Best effort for oblique
        si_in_vert = dir_cos(3, vert_img_dim_original);
        if si_in_vert > 0
            plotMap_display = flip(plotMap_display, 1);
            display_info.flipped_dims = [display_info.flipped_dims, 1];
        end
        
        lr_in_horiz = dir_cos(1, horiz_img_dim_original);
        if strcmp(convention, 'radiological')
            if lr_in_horiz > 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        else
            if lr_in_horiz < 0
                plotMap_display = flip(plotMap_display, 2);
                display_info.flipped_dims = [display_info.flipped_dims, 2];
            end
        end
end

%% Calculate montage layout
num_slices = idx_end - idx_start + 1;
if num_slices < 5
    slices_per_row = num_slices;
else
    slices_per_row = 5;
end

%% Get tile dimensions
tile_height = size(plotMap_display, 1);
tile_width = size(plotMap_display, 2);

%% Display montage
max_val = prctile(plotMap_display(:), 95);

montage(squeeze(plotMap_display(:, :, idx_start:idx_end)), ...
    'ThumbnailSize', [tile_height, tile_width], ...
    'Size', [ceil(num_slices/slices_per_row), slices_per_row], ...
    'DisplayRange', [0 max_val]);

axis image;
hold on;

%% Setup colormap and colorbar
if viridis_map
    colormap viridis
else
    colormap gray
end

clim([0 max_val]);
axis tight;
axis off;

%% Add target label
text(1, 5, target,...
    'Rotation', 0, 'Color', 'w', 'FontSize', 15,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

%% Add orientation labels
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

%% Setup axes and colorbar
set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.01, 0.05, 0.98, 0.9]);

if cbar
    cbar_handle = colorbar;
    set(cbar_handle, 'Color', [1 1 1]);
    cbar_handle.Label.String = target;
    cbar_handle.Label.Color = [1 1 1];
end

hold off;

end