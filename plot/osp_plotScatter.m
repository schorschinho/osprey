function [out_scat] = osp_plotScatter(MRSCont,model,quant,metab,corrData,corrDataName,VoxelIndex)
%% [out_scat] = osp_plotScatter(MRSCont,model,quant,metab,corrData,corrDataName)
% Creates correlation figure of  the chosen quantifcation and metabolite
% one figure contains a correlation analysis with subgroup correlations.
% If no groups are defined the distribution of the whole dataset will be shown. If no
% correlation measures are passed no correlation will be plotted.
%
%   USAGE:
%       [out_scat] = osp_plotScatter(MRSCont,which,metab,corrData,corrDataName,tit,GUI)
%
%   OUTPUTS:
%       [out_scat] = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       model    = Fitting Style
%       quant    = Quantification
%                   OPTIONS:    'tCr' (default)
%                               'rawWaterScaled'
%       quant    = Quantification
%       corrData = Data for correlation analysis
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2019-11-14)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is uses a modified version of the regression_line_ci
%       package to plot confidence intervals
%       Boris Gutman
%       https://de.mathworks.com/matlabcentral/fileexchange/39339-linear-regression-confidence-interval
%       The colorbrewer package is included for nicer colors
%       Charles
%       https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
%
%   HISTORY:
%       2019-11-12: First version of the code.
%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin<7
    VoxelIndex = 1;
    if nargin<6
        corrDataName = 'correlation measure';
        if nargin<4
            error('No correlation data passed. Please add correlation data to the MRSCont.');
            if nargin<3
                metab = 'GABA';
                if nargin<2
                    which = 'tCr';
                    if nargin<1
                        error('ERROR: no input Osprey container specified.  Aborting!!');
                    end
                end
            end
        end
    end
end

% Check that OspreyOverview has been run before
if ~MRSCont.flags.didOverview
    error('Trying to create overview plots, but no overview data has been created. Run OspreyOverview first.')
end
[~] = osp_Toolbox_Check ('OspreyProcess',MRSCont.flags.isGUI);

%%% 2. CREATE COLORMAP %%%
[cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
temp = cb(3,:);
cb(3,:) = cb(4,:);
cb(4,:) = temp;

%%% 3. EXTRACT METABOLITE CONCENTRATIONS%%%
if strcmp(quant,'AlphaCorrWaterScaled') || strcmp(quant,'AlphaCorrWaterScaledGroupNormed')
    idx_1  = 1;
    ConcData = MRSCont.quantify.tables.(model).(quant).(['Voxel_' num2str(VoxelIndex)]) {:,idx_1};
else
    idx_1  = find(strcmp(MRSCont.quantify.metabs.(model),metab));
    ConcData = MRSCont.quantify.tables.(model).(quant).(['Voxel_' num2str(VoxelIndex)]) {:,idx_1};
end
metabFlag = 0;
if isempty(ConcData)
    dim = size(ConcData);
    ConcData = zeros(dim(1),1);
end

if ischar(corrData)
    metabFlag = 1;
    idx_1  = find(strcmp(MRSCont.quantify.metabs.(model),corrData));
    corrData = MRSCont.quantify.tables.(model).(quant).(['Voxel_' num2str(VoxelIndex)]) {:,idx_1};
end

if strcmp(quant, 'tCr')
    ylab = [metab ' / tCr'];
    if metabFlag
        xlab = [corrDataName ' / tCr'];
    else
        xlab = corrDataName;
    end
end
if strcmp(quant, 'rawWaterScaled')
    ylab = [metab ' rawWaterScaled  (i.u.)'];
    if metabFlag
        xlab = [corrDataName ' rawWaterScaled  (i.u.)'];
    else
        xlab = corrDataName;
    end
end
if strcmp(quant, 'CSFWaterScaled')
    ylab = [metab ' CSFWaterScaled  (i.u.)'];
    if metabFlag
        xlab = [corrDataName ' CSFWaterScaled  (i.u.)'];
    else
        xlab = corrDataName;
    end
end
if strcmp(quant, 'TissCorrWaterScaled')
    ylab = [metab ' TissCorrWaterScaled  (i.u.)'];
    if metabFlag
        xlab = [corrDataName ' TissCorrWaterScaled  (i.u.)'];
    else
        xlab = corrDataName;
    end
end
if strcmp(quant, 'AlphaCorrWaterScaled')
    ylab = [metab ' AlphaCorrWaterScaled  (i.u.)'];
    if metabFlag
        xlab = [corrDataName ' AlphaCorrWaterScaled  (i.u.)'];
    else
        xlab = corrDataName;
    end
end
if strcmp(quant, 'AlphaCorrWaterScaledGroupNormed')
    ylab = [metab ' AlphaCorrWaterScaledGroupNormed  (i.u.)'];
    if metabFlag
        xlab = [corrDataName ' AlphaCorrWaterScaledGroupNormed  (i.u.)'];
    else
        xlab = corrDataName;
    end
end
%%% 4. CREATE CORRELATION PLOT %%%
% Generate a new figure and keep the handle memorized
out_scat = figure('Color', 'w');
% Scatter plot with separate groups and correlatios and create legend
for g = 1 : MRSCont.overview.NoGroups
    x_tmp = corrData(MRSCont.overview.groups == g);
    y_tmp = ConcData(MRSCont.overview.groups == g);
    scatter (x_tmp,y_tmp,'SizeData', 10, 'MarkerFaceColor', cb(g,:), 'MarkerEdgeColor', 'none');
    hold on
end
legend(MRSCont.overview.groupNames);
legend('boxoff');
legend('AutoUpdate','off','Location','north','Orientation','horizontal');
str = cell(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    x_tmp = corrData(MRSCont.overview.groups == g);
    y_tmp = ConcData(MRSCont.overview.groups == g);
    X = [ones(length(x_tmp),1) x_tmp];
    b = X\y_tmp;
    hold on
    [~, ~,~] = regression_line_ci(0.05,b,x_tmp,y_tmp,cb(g,:));
    [r,p] = corrcoef(x_tmp,y_tmp);
    if length(r) > 1
        str{g} = ['r = ' num2str(r(1,2),'%3.3g') '\newline' 'p = ' num2str(p(1,2),'%3.3g')];
    else
        str{g} = ['r = nan \newline' 'p = nan']; 
    end
end

% Scatter plot and correlation including all datapoints
if g ~= 1
    X = [ones(length(corrData),1) corrData];
    b = X\ConcData;
    [r,p] = corrcoef(corrData,ConcData);
    str{MRSCont.overview.NoGroups+1} = ['r = ' num2str(r(1,2),'%3.3g') '\newline' 'p = ' num2str(p(1,2),'%3.3g')];
    [~, ~,~] = regression_line_ci(0.05,b,corrData,ConcData,[0,0,0]);
end

%Adjusting size of the axis and setting up text for r and p values
cXlim = get(gca,'XLim');
cYlim = get(gca,'YLim');
magy = abs(cYlim(2) - cYlim(1));
set(gca, 'YLim', [cYlim(1)-(0.1* magy)  cYlim(2)]);
magx = abs(cXlim(2) - cXlim(1));
spacing = magx/length(str);
for an = 1 : length(str)
    if an == MRSCont.overview.NoGroups+1
                text(cXlim(1) + magx*0.01 + (spacing*(an-1)) , cYlim(1),str{an},'Color',[0,0,0]);
    else
        text(cXlim(1) + magx*0.01 + (spacing*(an-1)) , cYlim(1),str{an},'Color',cb(an,:));
    end
end

% Black axes, white background
if ~MRSCont.flags.isGUI
    set(gca, 'YColor', 'w');
    title([corrDataName ' vs ' metab],'FontSize',16);
else
    set(gca, 'YColor', MRSCont.colormap.Foreground);
    set(gca, 'XColor', MRSCont.colormap.Foreground);
    title([corrDataName ' vs ' metab],'FontSize',16, 'Color', MRSCont.colormap.Foreground);
end
box off;
xlabel(xlab,'FontSize',16);
ylabel(ylab,'FontSize',16);

%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out_scat, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end
