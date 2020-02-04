function out = osp_plotMeanSpec(MRSCont, which,GUI, g, shift,group, xlab, ylab, figTitle)
%% out = osp_plotMeanSpec(MRSCont, which,g, shift, GUI, xlab, ylab, figTitle)
%   Creates a figure mean and standard deviation of the spectra. If the
%   chosen spectra was fitted the mean fit, baseline and residue are shown.
%
%   USAGE:
%       out = osp_plotMeanSpec(MRSCont, which,GUI, g, shift,group, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'A' (default)
%                               'B'
%                               'C'
%                               'D'
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%       GUI       = flag to decide whether plot is used in GUI
%       g         = group index
%       shift     = shift between groups
%       group     = control flag to decide whether different groups are plotted (removes residual and baseline)  
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
if nargin<8
ylab='';
    if nargin<7
        xlab='Frequency (ppm)';
        if nargin<6
            group = 0;
            if nargin<5
                shift = 0.1;
                if nargin<4
                    g = 1;    
                    if nargin<3
                        GUI = 0;    
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
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract normalized spectra and fits
shift = shift * (g-1);
if MRSCont.flags.isUnEdited
    switch which
        case 'A'
            fit = 'off';
            fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
            fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
            data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
            data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
            baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
            baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
            residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
            residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);
            ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
        case {'ref','w'}
            fit = which;
            fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
            fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]);
            data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
            data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
            residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
            residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);             
            ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
    end
end
if MRSCont.flags.isMEGA
    switch which
        case 'A'
            name = 'off';
            if strcmp(fitStyle,'Concatenated')
                data_mean = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' which]);
                data_sd = MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' which]);
                ppm = MRSCont.overview.(['ppm_data_' which]);
            else
                fit = 'off';
                fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
                fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
                data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
                data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
                baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
                baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
                residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
                residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);
                ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
            end            
        case 'B'
            name = 'on';
            data_mean = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' which]);
            data_sd = MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' which]);  
            ppm = MRSCont.overview.(['ppm_data_' which]);
        case 'diff1'
            if ~strcmp(fitStyle,'Concatenated')
                fit = which;
            else
                fit = 'conc';
            end
            fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
            fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
            data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
            data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
            baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
            baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
            residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
            residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);
            ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
        case {'sum'}
            if ~strcmp(fitStyle,'Concatenated')
                fit = which;
                data_mean = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' which]);
                data_sd = MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' which]);  
                ppm = MRSCont.overview.(['ppm_data_' which]);
            else
                fit = 'conc';
                fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
                fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
                data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
                data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
                baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
                baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
                residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
                residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);
                ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
            end             
        case {'ref','w'}
            fit = which;
            fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
            fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]);
            data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
            data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
            residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
            residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);            
            ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
    end
end
if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
    switch which
        case {'A','B','C','D'}
            data_mean = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' which]);
            data_sd = MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' which]);
            ppm = MRSCont.overview.(['ppm_data_' which]);
        case {'diff1','diff2','sum'}
            if ~strcmp(fitStyle,'Concatenated')
                fit = which;
                fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
                fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
                data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
                data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
                baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
                baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
                residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
                residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]); 
                ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
            else
                fit = 'conc';
                fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
                fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]); 
                data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
                data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
                baseline_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' fit '_' which]);
                baseline_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' fit '_' which]);
                residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
                residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);
                ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
            end            
        case {'ref','w'}
            fit = which;
            fit_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' fit '_' which]);
            fit_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' fit '_' which]);
            data_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_data_' fit '_' which]);
            data_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_data_' fit '_' which]);
            residual_mean = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' fit '_' which]);
            residual_sd = MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' fit '_' which]);             
            ppm = MRSCont.overview.(['ppm_fit_' fit '_' which]);
    end
end

if ~exist('name', 'var')
    name = which;
end

if nargin<8    
    if (~strcmp(which,'w') && ~strcmp(which,'ref'))
        if exist('fit_mean', 'var')
            figTitle = ['mean data \pm SD & mean fit: ' which]; 
        else
            figTitle = ['mean \pm  SD: ' name];
        end
        ppmRange = MRSCont.opts.fit.range;
    else
       figTitle = ['mean data \pm SD & mean fit: ' which];
       ppmRange = MRSCont.opts.fit.rangeWater;
    end
end

%Calculate SD shadows
data_yu = data_mean + data_sd;
data_yl = data_mean - data_sd;


if exist('fit_mean', 'var')
    fit_yu = fit_mean + fit_sd;
    fit_yl = fit_mean - fit_sd;
end
if exist('baseline_mean', 'var')
    baseline_yu = baseline_mean + baseline_sd;
    baseline_yl = baseline_mean - baseline_sd;
end
if exist('residual_mean', 'var')
    residual_yu = residual_mean + residual_sd;
    residual_yl = residual_mean - residual_sd;
end
maxshift = max(data_yu);
maxshift_abs = max(abs(data_yu));
shift = maxshift_abs * shift;
%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
out = figure('Visible','on');
hold on

%%% 4. PLOT DATA, FIT, RESIDUAL, BASELINE %%%
% Determine a positive stagger to offset data, fit, residual, and 
% baseline from the individual metabolite contributions
if GUI
    if exist('residual_mean', 'var')
        if ~group
            fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
        end
    end
    fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data

    if exist('fit_mean', 'var')
        plot(ppm,fit_mean+shift ,'color', MRSCont.colormap.Accent, 'LineWidth', 1.5); %Fit
        if ~group
            plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))) ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);  %Residual
        end
    end

    if exist('baseline_mean', 'var')
        plot(ppm,baseline_mean+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1); %Baseline
    end

    plot(ppm,data_mean+shift ,'color',cb(g,:), 'LineWidth', 2); % Data 


    if exist('fit_mean', 'var')
        if ~group
            plot(ppm, (zeros(1,length(ppm))), 'Color',MRSCont.colormap.Foreground); % Zeroline
            text(ppm(1)-0.05, 0, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))), 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue 
            text(ppm(1)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
            text(ppm(1)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
        end
    end
else
    if exist('fit_mean', 'var')
        if ~group
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color','k', 'LineWidth', 1); % Maximum Data
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))), 'Color','k', 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), 'Color','k', 'LineWidth', 1); % Max Residue 
            text(ppm(1)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color','k'); %Zeroline Residual text
            text(ppm(1)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color','k'); %Max Residue text
        end
    end
    if exist('residual_mean', 'var')
        if ~group
            fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
        end
    end
    fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data 
    plot(ppm,data_mean+shift ,'color','k', 'LineWidth', 2); % Data 
    if exist('fit_mean', 'var')
        plot(ppm,fit_mean+shift ,'color', 'r', 'LineWidth', 1); %Fit
        if ~group
            plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))) ,'color', 'k', 'LineWidth', 1);  %Residual 
            plot(ppm, (zeros(1,length(ppm))), 'Color','k'); % Zeroline
            text(ppm(1)-0.05, 0, '0', 'FontSize', 10,'Color','k'); %Zeroline text  
        end
    end    
    if exist('baseline_mean', 'var')
        if ~group
            plot(ppm,baseline_mean+shift ,'color', 'b', 'LineWidth', 1); %Baseline
        end
    end                                        
  
end

%%% 5. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)]);
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    if ~GUI
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


%%% 6. ADD OSPREY LOGO %%%
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end

   