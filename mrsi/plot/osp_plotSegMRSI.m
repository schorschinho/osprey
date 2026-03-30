function out = osp_plotSegMRSI(MRSCont,kk,idx_start,idx_end, plot_mask)
%% out = osp_plotSegMRSI(MRSCont,kk,idx_start,idx_end, plot_mask)
%   Creates a figure showing partial volume fractions of the MRSI voxels.
%   If plot_mask is true it will also show the binary masks for brain and
%   lipid voxels.
%
%   USAGE:
%       out = osp_plotSegMRSI(MRSCont,kk,idx_start,idx_end, plot_mask)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       idx_start = index MRSI slice to start (1 at bottom)
%       idx_end = index MRSI slice to end
%       plot_mask = plot binary masks for lipids and brain
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

%% Fall back to defaults if not provided
if nargin < 5
    plot_mask = 0;
    if nargin < 4
        idx_end = MRSCont.raw{1, 1}.nZvoxels;
        if nargin < 3
            idx_start = 1;
            if nargin < 2
               kk = 1; 
                    if nargin<1
                        error('ERROR: no input Osprey container specified.  Aborting!!');
                    end
            end
        end
    end
end


%%  Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot coregistration, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot segmentation, but coregistration has not been performed yet. Run OspreyCoreg first.')
end

if ~MRSCont.flags.didSeg
    error('Trying to plot segmentation, but segmentation has not been performed yet. Run OspreySeg first.')
end

%% Setup figure
out = figure;  
if ~plot_mask
    tiledlayout(4,1,'TileSpacing','compact')
else
    tiledlayout(6,1,'TileSpacing','compact')
end

set(out, 'Color', [0 0 0]);


if idx_end - idx_start + 1 < 5
    slices_per_row = idx_end - idx_start + 1;
else
    slices_per_row = 5;
end

fraction_names = {'fGM','fWM','fCSF','fLIP'};
for i = 1 : 4
    fraction_img = squeeze(MRSCont.seg.tissue.(fraction_names{i})(kk,:,:,:))*100;
    if (MRSCont.raw{1}.nZvoxels > 1)
        fraction_img = flip(fraction_img,3);
    end
    fraction_img = flip(fraction_img,1);
    nexttile
    montage(rot90(squeeze(fraction_img(:, :, idx_start:idx_end))),...
                 'ThumbnailSize', [size(fraction_img, 2), size(fraction_img, 1)],...
                 'Size', [NaN slices_per_row],'DisplayRange',[0 100]);

    text(10, 3, fraction_names{i},...
        'Rotation', 0,'Color','w','FontSize',15,...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    if i == 1
        cbar = colorbar;
        cbar.Label.String = 'tissue fraction (%)';
        cbar.Color = [1 1 1];
        cbar.Location = 'north';
    end
    colormap('gray');
end

set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.01, 0.05, 0.98, 0.9]); % Example: 10% margin on all sides
current_pos = cbar.Position;
cbar.Position = [current_pos(1), current_pos(2)*1.1, current_pos(3), current_pos(4)];

if plot_mask
    fraction_names = {'brain','lip'};
    for i = 1 : 2
        fraction_img = squeeze(MRSCont.seg.tissue.(fraction_names{i})(kk,:,:,:));
        if (MRSCont.raw{1}.nZvoxels > 1)
            fraction_img = flip(fraction_img,3);
        end
        fraction_img = flip(fraction_img,1);
        nexttile
        montage(rot90(squeeze(fraction_img(:, :, idx_start:idx_end))),...
                     'ThumbnailSize', [size(fraction_img, 2), size(fraction_img, 1)],...
                     'Size', [NaN slices_per_row],'DisplayRange',[0 1]);
    
        text(10, 3, [ fraction_names{i} ' mask'],...
            'Rotation', 0,'Color','w','FontSize',15,...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        colormap('gray');
    end
end




end

   