function [out_rain] = osp_plotRaincloud(MRSCont,model,quant,metab,tit,GMean,VoxelIndex)
%% [out_rain] = osp_plotRaincloud(MRSCont,metab,plots,tit)
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
%       model    = Fitting style 
%       quant    = Quantification
%                   OPTIONS:    'tCr'
%                               'rawWaterScaled'
%       metab    = metabolite for analysis
%       tit      = Title of the raincloud plot
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2019-11-14)
%       hzoelln2@jhmi.edu
%
%   CREDITS:    
%       This code uses a modified version of the raincloud plots by Davide
%       Poggiali
%       https://github.com/RainCloudPlots/RainCloudPlots
%       DOI: 10.12688/wellcomeopenres.15191.1
%       The colorbrewer package is included for nicer colors
%       Charles
%       https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
%
%   HISTORY:
%       2019-11-12: First version of the code.

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin < 7
    VoxelIndex = 1;
    if nargin < 6
        GMean = 0;
        if nargin<5
            tit = '';
            if nargin<4
                metab = 'GABA';
                if nargin<3
                    quant = 'tCr';
                    if nargin<2
                        model = '';
                        if MRSCont.flags.isUnEdited
                            model = 'off';
                        end
                        if MRSCont.flags.isMEGA && strcmp(MRSCont.opts.fit.style,'Separate')
                                model = 'diff1';
                            else
                                model = 'conc';
                        end
                        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) && strcmp(MRSCont.opts.fit.style,'Separate')
                                model = 'diff1';
                            else
                                model = 'conc';
                        end
                        if strcmp(model,'')
                            error('ERROR: unable to retrieve default model, please specify.  Aborting!!');
                        end
                        if nargin<1
                            error('ERROR: no input Osprey container specified.  Aborting!!');
                        end
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

%%% 2. CREATE COLORMAP %%%
[cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
temp = cb(3,:);
cb(3,:) = cb(4,:);
cb(4,:) = temp;

%%% 3. EXTRACT METABOLITE CONCENTRATIONS OR QM%%%
if ~strcmp(quant,'Quality')
    if strcmp(quant,'AlphaCorrWaterScaled') || strcmp(quant,'AlphaCorrWaterScaledGroupNormed')
        idx_1  = 1;
        ConcData = MRSCont.quantify.tables.(model).(quant).(['Voxel_' num2str(VoxelIndex)]) {:,idx_1};  
    else
        idx_1  = find(strcmp(MRSCont.quantify.metabs.(model),metab));
        ConcData = MRSCont.quantify.tables.(model).(quant).(['Voxel_' num2str(VoxelIndex)]) {:,idx_1};  
    end
else
   quality = {'SNR','FWHM','freqShift'};
   idx_1  = find(strcmp(quality,metab));
   ConcData = MRSCont.QM.tables {:,idx_1};
   quality_Names = {'SNR','FWHM (ppm)','freqShift (Hz)'};
end

if isempty(ConcData)
    dim = size(ConcData);
    ConcData = zeros(dim(1),1);
end

if strcmp(quant, 'tCr')
    ylab = [metab ' / tCr'];
end
if strcmp(quant, 'rawWaterScaled')
    ylab = [metab ' rawWaterScaled  (i.u.)'];
end
if strcmp(quant, 'CSFWaterScaled')
    ylab = [metab ' CSFWaterScaled  (i.u.)'];
end
if strcmp(quant, 'TissCorrWaterScaled')
    ylab = [metab ' TissCorrWaterScaled  (i.u.)'];
end
if strcmp(quant, 'AlphaCorrWaterScaled')
    ylab = [metab ' AlphaCorrWaterScaled  (i.u.)'];
end
if strcmp(quant, 'AlphaCorrWaterScaledGroupNormed')
    ylab = [metab ' AlphaCorrWaterScaledGroupNormed  (i.u.)'];
end
if strcmp(quant, 'Quality')
    ylab = [quality_Names{idx_1}];
end
%%% 4. CREATE RAINCLOUD PLOT %%%
% Generate a new figure and keep the handle memorized
out_rain = figure('Color', 'w');
hold on
% Calculate ksdensities of all groups and keep the maximum for
% normalization
if ~GMean
    f = zeros(1,MRSCont.overview.NoGroups);
    for g = 1 : MRSCont.overview.NoGroups
        data{1,1} = ConcData(MRSCont.overview.groups == g);
        [f_tmp, ~, ~] = ksdensity(data{1,1}, 'bandwidth', []);
        f(g) =  max(f_tmp); 
    end
    maxim = max(f);
else
    f = zeros(1,1);
    data{1,1} = ConcData(:);
    [f_tmp, ~, ~] = ksdensity(data{1,1}, 'bandwidth', []);
    f(1) =  max(f_tmp); 
    maxim = max(f);    
end

% Calculate mean and SD for the plot and create legend
if ~GMean
    for g = 1 : MRSCont.overview.NoGroups
        data{1,1} = ConcData(MRSCont.overview.groups == g);
        % mean and SD
        meanv = mean(data{1,1});
        sdv = std(data{1,1});
        MaSD{g} = errorbar(meanv,.1*g,sdv,'horizontal');
        set(MaSD{g}, 'Color', [0 0 0], 'LineWidth', 1, 'Marker','diamond','MarkerEdgeColor',cb(g,:),'MarkerFaceColor',cb(g,:),  'CapSize', 8  )
        legend(MRSCont.overview.groupNames);
    end
else
    data{1,1} = ConcData(:);
    % mean and SD
    meanv = mean(data{1,1});
    sdv = std(data{1,1});
    MaSD{1} = errorbar(meanv,.1,sdv,'horizontal');
    set(MaSD{1}, 'Color', [0 0 0], 'LineWidth', 1, 'Marker','diamond','MarkerEdgeColor',cb(1,:),'MarkerFaceColor',cb(1,:),  'CapSize', 8  )
    legend({'Grand mean'});
end

legend('boxoff');
legend('AutoUpdate','off','Location','north','Orientation','horizontal');

% Create raincloud plots
if ~GMean
    for g = 1 : MRSCont.overview.NoGroups
        data{1,1} = ConcData(MRSCont.overview.groups == g);

        rain{g} = raincloud_plot(data{1,1}, 'box_on', 1, 'color', cb(g,:), 'alpha', 0.3,...
             'box_dodge', 1, 'box_dodge_amount', .15*g, 'dot_dodge_amount', .15*g,...
             'box_col_match', 1,'cloud_edge_col', 'none', 'normalize', maxim);
    end
else
    data{1,1} = ConcData(:);
    rain{1} = raincloud_plot(data{1,1}, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.3,...
         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
         'box_col_match', 1,'cloud_edge_col', 'none', 'normalize', maxim);
end

%Adjust y-axis
xlabel(ylab,'FontSize',16);
if ~GMean
    minYlim = rain{MRSCont.overview.NoGroups}{3}.Position(2);
else
     minYlim = rain{1}{3}.Position(2);
end
cYlim = get(gca,'YLim');
set(gca, 'YLim', [(minYlim - 0.05) 1.10]);

% Black axes, white background
if ~MRSCont.flags.isGUI
    set(gca, 'YColor', 'w');
    title([tit ': ' model ' ' metab],'FontSize',16);
else
    set(gca, 'YColor', MRSCont.colormap.Background);
    set(gca, 'XColor', MRSCont.colormap.Foreground);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})
    title([tit ': ' model ' ' metab],'FontSize',16, 'Color', MRSCont.colormap.Foreground);
end
box off

%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out_rain, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end