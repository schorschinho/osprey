function out = osp_plotOverviewMMSpec(MRSCont, which, g, shift, xlab, ylab, figTitle)
%% out = osp_plotOverviewSpec(MRSCont, which, g, shift, xlab, ylab, figTitle)
%   Creates a figure with all spectra or fits.
%
%   USAGE:
%       out = osp_plotOverviewSpec, which,GUI, g, shift, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       which    = String for the spectrum or fit to plot (optional)
%                   OPTIONS:    'A' (default)
%                               'B'
%                               'C'
%                               'D'
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%                               'Fit:off'
%                               'Fit:diff1'
%                               'Fit:diff2'
%                               'Fit:sum'
%                               'Fit:ref'
%                               'Fit:w'
%       g         = group index
%       shift     = shift between groups
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2019-10-02)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2020-01-24: First version of the code.

% Check that OspreyOverview has been run before
if ~MRSCont.flags.didOverview
    error('Trying to create overview plots, but no overview data has been created. Run OspreyOverview first.')
end

cb = cbrewer('qual', 'Dark2', 12, 'pchip');
temp = cb(3,:);
cb(3,:) = cb(4,:);
cb(4,:) = temp;

%%% 1. PARSE INPUT ARGUMENTS %%%
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
if nargin<7
ylab='';
    if nargin<6
        xlab='Frequency (ppm)';
        if nargin<5
            group = 0;
            if nargin<4
                shift = 0.1;
                if nargin<3
                    g = 1;      
                    if nargin < 2
                        which = 'A';
                        if nargin<1
                            error('ERROR: no input Osprey container specified.  Aborting!!');
                        end
                    end
                end
            end
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract normalized spectra and fits
shift = shift * (g-1);
%Is spectrum?
if (strcmp(which,'A') || strcmp(which,'B') || strcmp(which,'C') || strcmp(which,'D') || strcmp(which,'diff1') || strcmp(which,'diff2') || strcmp(which,'sum') || strcmp(which,'ref') || strcmp(which,'w'))
    data = MRSCont.overview.sort_data.(['g_' num2str(g)]).(which);
    if nargin<8    
        if (~strcmp(which,'w') && ~strcmp(which,'ref'))
            figTitle = ['Individual specs: ' which]; 
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['Individual specs: ' which]; 
           ppmRange = MRSCont.opts.fit.rangeWater;
        end
    end
else % Is fit?
    fitwhich = which(5:end);
    if MRSCont.flags.isUnEdited %Is UnEdited
        switch fitwhich
            case 'off'
                fit = 'A';
                data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fitwhich '_' fit]);
            case 'ref'
                fit = fitwhich;
                data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fitwhich '_' fit]);
        end
    end
    if MRSCont.flags.isMEGA %Is MEGA
        switch fitwhich
            case {'diff1','sum'}
                if strcmp(fitStyle,'Concatenated') %Is Concatenated?
                    fit = 'conc';
                    data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fit '_' fitwhich]);
                else
                    fit = 'diff1';
                    data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fitwhich '_' fit]);            
                end   
            case {'off'}
                    fit = 'A';
                    data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fitwhich '_' fit]);                         
            case {'ref','w'}
                fit = fitwhich;
                data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fit '_' fitwhich]);
        end
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is Multiplexed
        switch fitwhich
            case {'diff1','diff2','sum'}
                if strcmp(fitStyle,'Concatenated') %Is Concatenated?
                    fit = 'conc';
                    data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fit '_' fitwhich]);
                else
                    fit = fitwhich;
                    data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fit '_' fitwhich]);            
                end           
            case {'ref','w'}
                fit = fitwhich;
                data = MRSCont.overview.sort_fit.(['g_' num2str(g)]).([fitwhich '_' fit]);
        end
    end
    
    if nargin<8    
        if (~strcmp(fitwhich,'w') && ~strcmp(fitwhich,'ref'))
            figTitle = ['Individual fits: ' fit ' ' fitwhich]; 
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['Individual fits: ' fit ' ' fitwhich]; 
           ppmRange = MRSCont.opts.fit.rangeWater;
        end
    end
end
            


%%% 3. PLOT DATA %%%
if length(which)>4
    if ~strcmp(which(1:4),'Fit:')
        maxshift_abs = max(abs(data{1}.specs));
        shift = maxshift_abs * shift;
        out = plot(data{1}.ppm,data{1}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        hold on;
        for kk = 2 : length(data) 
            plot(data{kk}.ppm,data{kk}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        end
    else
        maxshift_abs = max(abs(data{1}.fit));
        shift = maxshift_abs * shift;
        out = plot(data{1}.ppm,data{1}.fit+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
        hold on;
        for kk = 2 : length(data) 
            plot(data{kk}.ppm,data{kk}.fit+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
        end
    end
else
        maxshift_abs = max(abs(data{1}.specs));
        shift = maxshift_abs * shift;
        out = plot(data{1}.ppm,data{1}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        hold on;
        for kk = 2 : length(data) 
            plot(data{kk}.ppm,data{kk}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        end
end
%%% 4. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)], 'XMinorTick', 'On');
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    if ~MRSCont.flags.isGUI
        set(gca, 'YColor', 'w');
        % Black axes, white background
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        set(gcf, 'Color', 'w');
        title(figTitle, 'Interpreter', 'Tex');
    else
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        title(figTitle, 'Interpreter', 'Tex', 'Color', MRSCont.colormap.Foreground);
    end
else
    set(gca, 'YColor', 'k');
end

box off;
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out.Parent.Parent, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end

   