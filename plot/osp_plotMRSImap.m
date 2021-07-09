function out = osp_plotMRSImap(MRSCont, kk, nominator_spec,denominator_spec, nominator, denominator, upsample,figTitle)
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
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
if nargin < 7
    upsample = 1;
    if nargin<6
        denominator = [];
        if nargin<5
            nominator = {'Cr', 'PCr'}; 
            if nargin<4
                denominator_spec = 'none'; 
                if nargin < 3
                    nominator_spec = 'off';
                    if nargin < 2
                        kk = 1;
                        if nargin<1
                            error('ERROR: no input Osprey container specified.  Aborting!!');
                        end
                    end
                end
            end
        end
    end
end
if nargin<8    
    [~,filen,ext] = fileparts(MRSCont.files{kk});
    nominator_name = '';
    denominator_name = '';
    for t = 1 : length(nominator)
        nominator_name = [nominator_name nominator{t}];
        if t < length(nominator)
            nominator_name = [nominator_name '+'];
        end
    end
    if ~isempty(denominator)
        for t = 1 : length(denominator)
            denominator_name = [denominator_name denominator{t}];
            if t < length(denominator)
            denominator_name = [denominator_name '+'];
            end
        end
        figTitle = sprintf([filen ext '\n' fitMethod ' ' fitStyle ': ' nominator_name '/' denominator_name  ]); 
    else
        figTitle = sprintf([filen ext '\n' fitMethod ' ' fitStyle ': ' nominator_name ' raw amplitudes']); 
    end    
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
nominator_map = zeros(size(MRSCont.quantify.amplMets{kk}));
denominator_map = zeros(size(MRSCont.quantify.amplMets{kk}));
if  (MRSCont.flags.isMRSI == 1)
    for nom = 1 : length(nominator)
        nominator_map = nominator_map + MRSCont.quantify.amplMets{kk}.(nominator_spec).(nominator{nom});
    end
    if ~isempty(denominator) && ~isempty(denominator_spec)
        for denom = 1 : length(denominator)
            denominator_map = denominator_map + MRSCont.quantify.amplMets{kk}.(denominator_spec).(denominator{denom});
        end
        map = nominator_map ./ denominator_map;
        map(denominator_map == 0) = 0;
    else
        map = nominator_map;
    end    
end

map_mean = mean(mean(map));

%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
canvasSize  = get(0,'defaultfigureposition');
if ~MRSCont.flags.isGUI
    out = figure('Position', canvasSize);
else
    out = figure('Position', canvasSize,'Visible','off');
end
map = map';

sz_map = size(map);

mask = MRSCont.mask{kk};
map = mask .* map';
if upsample > 1
    map = imresize(map,sz_map .* upsample);
end





colormap = viridis(100);
heatmap(map,'Colormap',colormap);


% heatmap(map,'Colormap',gray);

caxis(out.Children,[0 2]);
colorbar off

%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes

set(gca, 'FontSize', 16);
ax = gca;
% 
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

if ~MRSCont.flags.isGUI
    % Black axes, white background
    title(figTitle);
else
    title(figTitle);
end

set(gcf, 'Color', MRSCont.colormap.Background);   

set(gcf,'Units','Normalized');
set(gcf,'Position',[1.1250 0.1411 0.66 1]);
%%% 8. ADD OSPREY LOGO AND TIGHTEN FIGURE %%%
% if ~MRSCont.flags.isGUI
%     [I, map] = imread('osprey.gif','gif');
%     axes(out, 'Position', [0, 0.9, 0.1, 0.1*11.63/14.22]);
%     imshow(I, map);
%     axis off;
% end

end

   