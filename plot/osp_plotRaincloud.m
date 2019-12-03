function [out_rain] = osp_plotRaincloud(MRSCont,which,metab,tit,GUI)
%% [out_rain] = osp_plotRaincloud(MRSCont,which,metab,plots,tit,GUI)
% Creates raincloud plot from the quantification tables and the chosen quantifcation and metabolite
% The figure contains raincloud plots with boxplots, as well as, mean
% and sd of the indicated groups in the overview struct. If no groups are
% defined the distribution of the whole dataset will be shown.
%
%   USAGE:
%       [out_rain] = osp_plotRaincloud(MRSCont,which,metab,plots,corrData,corrDataName,tit,GUI)
%
%   OUTPUTS:
%       [out_rain] = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       which    = Quantification
%                   OPTIONS:    'tCr'
%                               'rawWaterScaled'
%       metab    = metabolite for analysis
%       tit      = Title of the raincloud plot
%       GUI      = flag if fiure is used in GUI
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2019-11-14)
%       hzoelln2@jhmi.edu
%
%   CREDITS:    
%       This code uses a modified version of the raincloud plots by David
%       Poggiali
%       https://github.com/RainCloudPlots/RainCloudPlots
%       DOI: 10.12688/wellcomeopenres.15191.1
%       The colorbrewer package is included for nicer colors
%       Charles
%       https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
%
%   HISTORY:
%       2019-11-12: First version of the code.


% Check that OspreyOverview has been run before
if ~MRSCont.flags.didOverview
    error('Trying to create overview plots, but no overview data has been created. Run OspreyOverview first.')
end

%%% 1. CREATE COLORMAP %%%
[cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
temp = cb(3,:);
cb(3,:) = cb(4,:);
cb(4,:) = temp;

%%% 2. EXTRACT METABOLITE CONCENTRATIONS%%%
idx_1  = find(strcmp(MRSCont.quantify.metabs,metab));
ConcData = MRSCont.quantify.tables.(which) {:,idx_1};      

if strcmp(which, 'tCr')
    ylab = [metab ' / tCr'];
end
if strcmp(which, 'rawWaterScaled')
    ylab = [metab ' rawWaterScaled  (i.u.)'];
end
if strcmp(which, 'CSFWaterScaled')
    ylab = [metab ' CSFWaterScaled  (i.u.)'];
end
if strcmp(which, 'TissCorrWaterScaled')
    ylab = [metab 'TissCorrWaterScaled  (i.u.)'];
end
%%% 3. CREATE RAINCLOUD PLOT %%%
% Generate a new figure and keep the handle memorized
out_rain = figure('Color', 'w');
hold on
% Calculate ksdensities of all groups and keep the maximum for
% normalization
f = zeros(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    data{1,1} = ConcData(MRSCont.overview.groups == g);
    [f_tmp, ~, ~] = ksdensity(data{1,1}, 'bandwidth', []);
    f(g) =  max(f_tmp); 
end
maxim = max(f);

% Calculate mean and SD for the plot and create legend
for g = 1 : MRSCont.overview.NoGroups
    data{1,1} = ConcData(MRSCont.overview.groups == g);
    % mean and SD
    meanv = mean(data{1,1});
    sdv = std(data{1,1});
    MaSD{g} = errorbar(meanv,.1*g,sdv,'horizontal');
    set(MaSD{g}, 'Color', [0 0 0], 'LineWidth', 1, 'Marker','diamond','MarkerEdgeColor',cb(g,:),'MarkerFaceColor',cb(g,:),  'CapSize', 8  )
end
legend(MRSCont.overview.groupNames);
legend('boxoff');
legend('AutoUpdate','off','Location','north','Orientation','horizontal');

% Create raincloud plots
for g = 1 : MRSCont.overview.NoGroups
    data{1,1} = ConcData(MRSCont.overview.groups == g);

    rain{g} = raincloud_plot(data{1,1}, 'box_on', 1, 'color', cb(g,:), 'alpha', 0.3,...
         'box_dodge', 1, 'box_dodge_amount', .15*g, 'dot_dodge_amount', .15*g,...
         'box_col_match', 1,'cloud_edge_col', 'none', 'normalize', maxim);
end

%Adjust y-axis
xlabel(ylab,'FontSize',16);
minYlim = rain{MRSCont.overview.NoGroups}{3}.Position(2);
cYlim = get(gca,'YLim');
set(gca, 'YLim', [(minYlim - 0.05) 1.10]);

% Black axes, white background
if ~GUI
    set(gca, 'YColor', 'w');
    title([tit ': ' metab],'FontSize',16);
else
    set(gca, 'YColor', MRSCont.colormap.Background);
    set(gca, 'XColor', MRSCont.colormap.Foreground);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})
    title([tit ': ' metab],'FontSize',16, 'Color', MRSCont.colormap.Foreground);
end
box off

%%% 4. ADD OSPREY LOGO %%%
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end