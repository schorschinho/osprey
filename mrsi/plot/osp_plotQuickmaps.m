function out = osp_plotQuickmaps(MRSCont, spec, target,idx_start,idx_end,viridis_map,interpolated,brainMask)
%% out = oosp_plotQuickmaps(MRSCont, spec, target,idx_start,idx_end,viridis_map,interpolated,brainMask)
%   Creates a figure showing the quick integration maps of an MRSI scan
%
%   USAGE:
%       out = osp_plotQuickmaps(MRSCont, spec, target,idx_start,idx_end,viridis_map,interpolated,brainMask)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   ARGUMENTS:
%       MRSCont  = Osprey data container.
%       spec     = Target spectrum
%       target   = Target metabolite ('tNAA')
%       idx_start = index MRSI slice to start (1 at bottom)
%       idx_end = index MRSI slice to end
%       viridis_map = use viridis map
%       interpolated = interpolation factor
%       brainMask = apply brain mask first
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

% Check that OspreyCoreg has been run before
if ~MRSCont.flags.didLoadData && (strcmp(spec,'raw') || strcmp(spec,'raw_w'))
    error('Trying to plot quick maps, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didProcess && (strcmp(spec,'A'))
    error('Trying to plot quick maps, but data has not been processed yet. Run OspreyProcess first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
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
                                    if nargin<1
                                        error('ERROR: no input Osprey container specified.  Aborting!!');
                                    end
                            end
                        end
                end
            end
        end
    end
end
if interpolated % This does not work together
    brainMask = 0;
end
if brainMask &&  ~MRSCont.flags.didSeg
    msg('Trying to apply a brain mask, but data has not been segmented yet. Run OspreySeg first.')
    brainMask = 0;
end


if ~MRSCont.flags.isGUI
    out = figure;   
else
    out = figure('Visible','off');
end

if ~viridis_map
    set(out, 'Color', [0 0 0]); 
else
    vir = viridis;
    set(out, 'Color', vir(1,:)); 
end

if interpolated
    plotMap = MRSCont.quickMapsInt.(spec).(target);
else
    plotMap = MRSCont.quickMaps.(spec).(target);
    switch brainMask
        case 1
            mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
            if (MRSCont.raw{1}.nZvoxels > 1)
                mask = flip(mask,3);
            end
            mask = flip(mask,1);
            plotMap = plotMap .* mask;
        case 2
            mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
            if (MRSCont.raw{1}.nZvoxels > 1)
                mask = flip(mask,3);
            end
            mask = flip(mask,1);
            plotMap = mask; % For debugging it might be useful to plot the brain mask only
        otherwise
    end
end

% For NIfTI convention we don't flip
% if (MRSCont.raw{1}.nZvoxels > 1)
%     plotMap = flip(plotMap,3);
% end

% plotMap = flip(plotMap,1);

map_cat = rot90(squeeze(plotMap(:, :, idx_start:idx_end)));
map_cat_reshape = reshape(map_cat,[size(map_cat,1) size(map_cat,2) * size(map_cat,3)]);
max_val = prctile(plotMap(:), 95);
imagesc(map_cat_reshape);
axis image
if viridis_map
    colormap viridis
else
    colormap gray
end

clim([0 max_val]) 
axis tight;
axis off;
hold on

text(1, 5, target,...
    'Rotation', 0,'Color','w','FontSize',15,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% set(out,'Units','normalized');
% pos = get(out, 'Position');
% 
set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.01, 0.05, 0.98, 0.9]); % Example: 10% margin on all sides
% set(out, 'Position', [pos(1) pos(2) (pos(4)/(map_cat_reshape(1)/map_cat_reshape(2))) pos(4)])
colorbar

end

   