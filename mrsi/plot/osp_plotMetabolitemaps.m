%% out = osp_plotMetabolitemaps(MRSCont, quantification, metabolite, idx_start, idx_end, interpolated, viridis_map, percentile_clean, clip, convention, cbar)
%   This function creates a figure showing metabolite maps from an MRSI
%   scan as a montage display.
%
%   The function handles axial, sagittal, coronal, and oblique acquisitions
%   by automatically detecting the image orientation and applying appropriate
%   transformations for correct anatomical display.
%
%   Metabolite maps can be displayed with optional interpolation, percentile
%   cleanup, and various colormap options including special QC visualization.
%
%
%   INPUTS:
%       MRSCont          = Osprey MRS data container.
%       quantification   = Quantification type string:
%                          'amplitudes', 'tCr', 'rawWaterScaled', 'CSFWaterScaled',
%                          'TissueCorr', 'AlphaCorr', 'water', 'GlobalQC', etc.
%       metabolite       = Target metabolite string ('tNAA', 'mI', 'Glx', etc.)
%                          Not used for 'water' or 'GlobalQC' quantification.
%       idx_start        = Starting MRSI slice index (optional, default: 1)
%       idx_end          = Ending MRSI slice index (optional, default: nZvoxels)
%       interpolated     = Interpolation factor (optional, default: 2)
%                          1 = no interpolation, >1 = interpolation scale factor
%       viridis_map      = Use viridis colormap (1) or gray colormap (0)
%                          (optional, default: 1)
%       percentile_clean = Apply percentile cleanup to remove outliers
%                          (0 = off, >0 = apply 97th percentile threshold)
%                          (optional, default: 0)
%       clip            = clip colormap (1 = clip, 0 = remove)
%       convention       = Display convention (optional, default: 'neurological'):
%                          'radiological' - Left on right side of image
%                          'neurological' - Left on left side of image
%       cbar             = Display colorbar (0 or 1) (optional, default: 1)
%
%   OUTPUTS:
%       out              = MATLAB figure handle.
%
%   PREREQUISITES:
%       - OspreyFit must have been run (MRSCont.flags.didFit = true)
%       - OspreyCoreg should have been run for orientation detection
%
%   NOTES:
%       - For QC maps, a special categorical colormap is used
%       - Complex values are automatically converted to real
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is part of the Osprey MRS analysis toolbox.
%       https://github.com/schorschinho/osprey
%
%   HISTORY:
%       2024-01-06: Added orientation labels and convention options.
%       2025-08-04: First version of the code.
%
%   SEE ALSO:
%       osp_plotMetabolitemapsOverlay, osp_plotQuickmaps,
%       osp_analyze_orientation, osp_add_orientation_labels_to_montage
%%

function out = osp_plotMetabolitemaps(MRSCont, quantification, metabolite, idx_start, idx_end, interpolated, viridis_map, percentile_clean,clip, convention, cbar)

%% Validate prerequisites
if ~MRSCont.flags.didFit
    error('Trying to plot metabolite maps, but data has not been modelled yet. Run OspreyFit first.')
end

%% Fall back to defaults if not provided
if nargin < 11
    cbar = 1;
    if nargin < 10
        convention = 'neurological';
      if nargin < 9  
        clip = 1;
        if nargin < 8
            percentile_clean = 0;
            if nargin < 7
                viridis_map = 1;
                if nargin < 6
                    interpolated = 2;
                    if nargin < 5
                        idx_end = MRSCont.raw{1, 1}.nZvoxels;
                        if nargin < 4
                            idx_start = 1;
                            if nargin < 3
                                metabolite = 'mI';
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

%% Validate convention input
if ~ismember(lower(convention), {'radiological', 'neurological'})
    warning('Invalid convention "%s". Using "neurological" as default.', convention);
    convention = 'neurological';
end
convention = lower(convention);

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
if ~strcmp(quantification,'water') && ~strcmp(quantification,'GlobalQC')
    plotMap = MRSCont.quantify.(quantification).(metabolite);
else
    plotMap = MRSCont.quantify.(quantification);
end

if ~isreal(plotMap)
    plotMap = real(plotMap);
end

%% Apply interpolation
if interpolated > 1
    plotMap = imresize3(plotMap, 'Scale', [interpolated interpolated 1], 'Method', 'cubic');
end

%% Get slices as requested (for percentile calculation before orientation transforms)
plotMapTemp = squeeze(plotMap(:, :, idx_start:idx_end));

%% Apply percentile filter if requested
if percentile_clean == 1
    max_val_clean = prctile(plotMapTemp(:), 97);
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

%% Prepare map for display with correct orientation
% Create a pseudo display_info structure for the map
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

%% Calculate display range
if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
    if percentile_clean == 1
            max_val = prctile(plotMapTemp, 97);
    elseif percentile_clean > 1
        max_val = percentile_clean;
    end
else
    max_val = max(plotMapTemp(:), [], 'all');
end

%% Display montage
montage(squeeze(plotMap_display(:, :, idx_start:idx_end)), ...
    'ThumbnailSize', [tile_height, tile_width], ...
    'Size', [ceil(num_slices/slices_per_row), slices_per_row], ...
    'DisplayRange', [0 max_val]);

axis image;
hold on;

%% Setup colormap
if viridis_map
    colormap viridis
else
    colormap gray
end
clim([0 max_val]); 

%% Handle QC-specific colormap
if contains(quantification,'QC') && ~contains(quantification,'QCfilt')  
    set(out, 'Color', [0 0 0]);
    if strcmp(quantification,'GlobalQC')
        mymap = [0 0 0
                1 0 0
                0 0 1
                0 1 0];
    else
        mymap = [0 0 0
                1 0 0           
                0 0 1
                1 0 1
                0 1 1
                0 1 0];
    end
    colormap(mymap);
end

axis tight;
axis off;

%% Add metabolite label
if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
    text(10, 3, metabolite,...
        'Rotation', 0, 'Color', 'w', 'FontSize', 15,...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

%% Add orientation labels
osp_add_orientation_labels_to_montage(idx_start, idx_end, slices_per_row, ...
    tile_width, tile_height, orientation_info, display_info);

%% Setup axes position
set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.01, 0.05, 0.98, 0.9]);

%% Add colorbar if requested
if cbar
    cbar_handle = colorbar;
    
    if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
        cbar_handle.Label.String = quantification;
    else
        cbar_handle.Label.String = 'Applied QC filter';
        if strcmp(quantification,'GlobalQC')
            cbar_handle.Ticks = [0 1 2 3];
            cbar_handle.TickLabels = {'brain mask', 'FWHM', 'SNR', 'QC passed'};
        else
            cbar_handle.Ticks = [0 1 2 3 4 5];
            cbar_handle.TickLabels = {'brain mask', 'FWHM', 'SNR', 'CRLB', 'percentile', 'QC passed'};
        end
    end
    cbar_handle.Color = [1 1 1];
end

hold off;

end