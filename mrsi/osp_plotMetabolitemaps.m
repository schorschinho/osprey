function out = osp_plotMetabolitemaps(MRSCont,quantification,metabolite,idx_start,idx_end,interpolated,viridis_map,percentile_clean)
%% out = osp_plotMetabolitemaps(MRSCont,quantification,metabolite,idx_start,idx_end,interpolated,viridis_map,percentile_clean)
%   Creates a figure showing the a metabolite map of an MRSI scan
%
%   USAGE:
%       out = osp_plotMetabolitemaps(MRSCont,quantification,metabolite,idx_start,idx_end,interpolated,viridis_map,percentile_clean)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       quantification = Which quantification to plot
%       metabolite   = Target metabolite ('tNAA')
%       idx_start = index MRSI slice to start (1 at bottom)
%       idx_end = index MRSI slice to end
%       interpolated = interpolation factor
%       viridis_map = use viridis map
%       percentile_clean = apply a percentile clean up, data outside is clipped
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot metabolite maps, but data has not been modelled yet. Run OspreyFit first.')
end



%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
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

if ~strcmp(quantification,'water') && ~strcmp(quantification,'GlobalQC')
    plotMap = MRSCont.quantify.(quantification).(metabolite);
else
    plotMap = MRSCont.quantify.(quantification);
end


% if (MRSCont.raw{1}.nZvoxels > 1)
%     plotMap = flip(plotMap,3);
% end


if interpolated > 1
    plotMap = imresize3(plotMap, 'Scale', [interpolated interpolated 1], 'Method' ,'cubic');
end

plotMapTemp = squeeze(plotMap(:, :, idx_start:idx_end));

if percentile_clean
    max_val = prctile(plotMapTemp(:), 97);
    plotMap(plotMap>max_val)=NaN;
end

map_cat = flip(rot90(flip(squeeze(plotMap(:, :, idx_start:idx_end)),3)),2);
map_cat_reshape = reshape(map_cat,[size(map_cat,1) size(map_cat,2) * size(map_cat,3)]);

if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
    max_val = prctile(plotMapTemp(:), 97);
else
    max_val = max(plotMapTemp(:),[], 'all');
end

imagesc(map_cat_reshape);
axis image

if viridis_map
    colormap viridis
else
    colormap gray
end

clim([0 max_val]); 

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
hold on

if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
    text(10, 3, metabolite,...
        'Rotation', 0,'Color','w','FontSize',15,...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','none');
end

set(gca, 'Units', 'normalized');
set(gca, 'Position', [0.01, 0.05, 0.98, 0.9]); % Example: 10% margin on all sides
cbar = colorbar;

if ~(contains(quantification,'QC') && ~contains(quantification,'QCfilt'))
    cbar.Label.String = quantification;
else
    cbar.Label.String = 'Applied QC filter';
    if strcmp(quantification,'GlobalQC')
        cbar.Ticks =  [0 1  2 3];
        cbar.TickLabels = {'brain mask', 'FWHM', 'SNR', 'QC passed'};
    else
        cbar.Ticks = [0  1 2 3 4 5];
        cbar.TickLabels = {'brain mask', 'FWHM', 'SNR', 'CRLB', 'percentile', 'QC passed'};
    end
end
cbar.Color = [1 1 1];
end

   