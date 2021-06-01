function out = osp_plotRawMRSIpos(MRSCont, kk, VoxelIndex)
%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%                                 'mm' re_mm
%
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didLoadData
    
error('Trying to visualize voxel position, but data has not been loaded yet. Run OspreyLoad first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin < 3
    VoxelIndex = [1 1];
    if nargin < 2
        kk = 1;
        if nargin<1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end

%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
%map = zeros(MRSCont.raw{kk}.nXvoxels, MRSCont.raw{kk}.nYvoxels);

% map = squeeze(abs(MRSCont.raw{kk}.fids(1,1,:,:)));
% map = map/(max(max(max(map))));
% map(map > 2 * map(1,1)) = 0.5;
% map(map < 0.5) = 0;

if MRSCont.flags.didSeg
        map = MRSCont.mask{kk};
        map = map(:,:,VoxelIndex(3));
        map = map'; 
        map(map > 0) = 0.5;
else if ~MRSCont.flags.didQuantify
    if ~MRSCont.flags.hasWater
        spec = op_freqrange(MRSCont.raw{kk},4.0,5.5);
        if MRSCont.flags.isMEGA
            spec = op_takesubspec(spec,1);
            spec                 = op_averaging(spec);            % Average
        end
        map = squeeze(sum(squeeze(abs(real(spec.specs))),1));
        map = map/(max(max(max(map))));
        map(map > 2 * map(1,1)) = 0.5;
        map(map < 0.5) = 0;
    else
        spec = op_freqrange(MRSCont.raw_w{kk},0,4.68*2);
%         if MRSCont.flags.isMEGA 
%             spec = op_takesubspec(spec,1);
%             spec                 = op_averaging(spec);            % Average
%         end
        map = squeeze(max(squeeze(abs(real(spec.specs)))));
        map = map/(max(max(max(map))));
        map(map > 0.05) = 0.5;
        map(map < 0.5) = 0;
    end
else    
     map = MRSCont.quantify.amplMets{kk}.off.NAA + MRSCont.quantify.amplMets{kk}.off.NAAG;
    map = map/(max(max(max(map))));
    map(map > 0.1) = 0.5;
    map(map < 0.5) = 0;
end
    
end


map(VoxelIndex(2),VoxelIndex(1)) = 1;
if size(map,3) <= 1
    map = map';
else
    map = squeeze(map(:,:,1))'
end

%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~MRSCont.flags.isGUI
    out = figure;
    set(gcf, 'Color', 'w');
else
    out = figure('Visible','off');
end

Background = [255/255 254/255 254/255];
Foreground = [11/255 71/255 111/255];
Accent = [254/255 186/255 47/255];
cmap = [Background; Foreground; Accent];
heatmap(map, 'CellLabelColor','none')
colorbar off
colormap(cmap);





%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes

set(gca, 'FontSize', 16);
ax = gca;

iniXLabel = get(ax,'XDisplayLabels');
iniYLabel = get(ax,'YDisplayLabels');

for l = 1 : length(iniXLabel)
        iniXLabel{l,1} = '';
end

for l = 1 : length(iniYLabel)
        iniYLabel{l,1} = '';
end
set(ax, 'XDisplayLabels',iniXLabel)
set(ax, 'YDisplayLabels',iniYLabel)
set(ax, 'FontSize',12)


set(gcf, 'Color', MRSCont.colormap.Background);   


%%% 8. ADD OSPREY LOGO AND TIGHTEN FIGURE %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.9, 0.1, 0.1*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end

   