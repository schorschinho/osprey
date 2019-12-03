function out = osp_plotProcess(MRSCont, kk, which, GUI, ppmmin, ppmmax)
%% out = osp_plotProcess(MRSCont, kk, which, GUI, ppmmin, ppmmax)
%   Creates a figure showing processed data stored in an Osprey data container,
%   ie in the raw fields. This function will display the *processed and
%   averaged* data, i.e. after spectral alignment, averaging, water removal,
%   and other processing steps carried out in OspreyProcess.
%
%   USAGE:
%       out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax, xlab, ylab)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to fit (optional)
%                   OPTIONS:    'A' (default)
%                               'B' (for MEGA, HERMES, HERCULES)
%                               'C' (for HERMES, HERCULES)
%                               'D' (for HERMES, HERCULES)
%                               'diff1' (for MEGA, HERMES, HERCULES)
%                               'diff2' (for HERMES, HERCULES)
%                               'sum' (for MEGA, HERMES, HERCULES)
%                               'ref'
%                               'w'
%       xlab     = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab     = label for the y-axis (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyProcess has been run before
if ~MRSCont.flags.didProcess
    error('Trying to plot processed data, but data has not been processed yet. Run OspreyProcess first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin<6
    switch which
        case {'mets'}
            ppmmax = 4.5;
        case {'ref', 'w'}
            ppmmax = 2*4.68;
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    if nargin<5
        switch which
            case {'mets'}
                ppmmin = 0.2;
            case {'ref', 'w'}
                ppmmin = 0;
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end
        if nargin < 4
            GUI = 0;
            if nargin < 3
                which = 'A';
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


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range
switch which
    case {'mets'}
        rawDataToPlot  = MRSCont.raw{kk};
        procDataToPlot = MRSCont.processed.A{kk};
        which = 'A';
    case 'ref'
        rawDataToPlot  = MRSCont.raw_ref{kk};
        procDataToPlot = MRSCont.processed.ref{kk};
    case 'w'
        rawDataToPlot  = MRSCont.raw_w{kk};
        procDataToPlot = MRSCont.processed.w{kk};
    otherwise
        error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
end


%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~GUI
    out = figure;
else
    out = figure('Visible','off');
end

nAvgs = rawDataToPlot.averages;
% Divide the figure into six tiles, create four axes
ax_raw      = subplot(3, 2, 1);
ax_aligned  = subplot(3, 2, 3);
ax_proc     = subplot(3, 2, 5);
ax_drift    = subplot(3, 2, 2);


%%% 4. PLOT RAW UNALIGNED %%%
% Add the data and plot
hold(ax_raw, 'on');    
% Loop over all averages
if ~GUI
    for rr = 1:nAvgs
        plot(ax_raw, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5);
    end
else
    for rr = 1:nAvgs
        plot(ax_raw, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5, 'Color', MRSCont.colormap.Foreground);
    end
end
plotRange = op_freqrange(rawDataToPlot, ppmmin, ppmmax);
yLims = [mean(min(real(plotRange.specs))) mean(max(real(plotRange.specs)))];
set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
y = ylim;
if ~(strcmp(which,'w') || strcmp(which,'ref'))
    if ~GUI
        plot(ax_raw, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_raw, [3.027 3.027], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_raw, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
    else
        plot(ax_raw, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', MRSCont.colormap.Foreground, 'LineWidth', 0.5);
        plot(ax_raw, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', MRSCont.colormap.Foreground, 'LineWidth', 0.5);
        plot(ax_raw, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', MRSCont.colormap.Foreground, 'LineWidth', 0.5);
    end
        
end
hold(ax_raw, 'off');
title(ax_raw, 'Pre-alignment', 'Color', MRSCont.colormap.Foreground);
xlabel(ax_raw, 'Frequency (ppm)', 'Color', MRSCont.colormap.Foreground)
if GUI
    set(ax_raw, 'YColor', MRSCont.colormap.Background);
    set(ax_raw,'YTickLabel',{})
    set(ax_raw,'YTick',{})
end


%%% 5. PLOT RAW ALIGNED %%%
% Apply stored corrections to calculate the spectra to display
applyDataToPlot = rawDataToPlot;
t = rawDataToPlot.t;
fs = procDataToPlot.specReg.fs;
phs = procDataToPlot.specReg.phs;
if isfield(MRSCont.QM.freqShift, which)
    refShift = -repmat(MRSCont.QM.freqShift.(which)(kk), size(fs));
    fs = fs + refShift;
end
for jj = 1:size(applyDataToPlot.fids,2)
    applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
end
applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);

hold(ax_aligned, 'on');    
% Loop over all averages
if ~GUI
    for rr = 1:nAvgs
        plot(ax_aligned, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5);
    end
else
    for rr = 1:nAvgs
        plot(ax_aligned, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5, 'Color', MRSCont.colormap.Foreground);
    end
end
set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
y = ylim;
if ~(strcmp(which,'w') || strcmp(which,'ref'))
    if ~GUI    
        plot(ax_aligned, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_aligned, [3.027 3.027], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_aligned, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
    else
        plot(ax_aligned, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_aligned, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_aligned, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);        
    end
end
hold(ax_aligned, 'off');
title(ax_aligned, 'Post-alignment', 'Color', MRSCont.colormap.Foreground);
xlabel(ax_aligned, 'Frequency (ppm)', 'Color', MRSCont.colormap.Foreground)
if GUI
    set(ax_aligned, 'YColor', MRSCont.colormap.Background);
    set(ax_aligned,'YTickLabel',{})
    set(ax_aligned,'YTick',{})
end


%%% 6. PLOT PROCESSED %%%
% Add the data and plot
hold(ax_proc, 'on');
if ~GUI 
    plot(ax_proc, procDataToPlot.ppm, procDataToPlot.specs/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'k', 'LineWidth', 1.5);
else
    plot(ax_proc, procDataToPlot.ppm, procDataToPlot.specs/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
end
y = [-0.2, 1.2];
set(ax_proc, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', y);
if ~(strcmp(which,'w') || strcmp(which,'ref'))
    if ~GUI    
        plot(ax_proc, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_proc, [3.027 3.027], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
        plot(ax_proc, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
    else
        plot(ax_proc, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_proc, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_proc, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground,  'LineWidth', 0.5);        
    end    
end
hold(ax_proc, 'off');
title(ax_proc, 'Aligned and averaged', 'Color', MRSCont.colormap.Foreground);
xlabel(ax_proc, 'Frequency (ppm)', 'Color', MRSCont.colormap.Foreground)
if GUI
    set(ax_proc, 'YColor', MRSCont.colormap.Background);
    set(ax_proc,'YTickLabel',{})
    set(ax_proc,'YTick',{})
end

%%% 7. GENERATE DRIFT PLOT %%%
if isfield(MRSCont.QM.drift, which)
    if length(MRSCont.QM.drift.(which){kk}) > 1
        crDriftPre = MRSCont.QM.drift.(which){kk};
        crDriftPost = op_measureDrift(applyDataToPlot);
        hold(ax_drift, 'on');
        if ~GUI
            plot(ax_drift, crDriftPre, 'o', 'Color', 'r');
            plot(ax_drift, crDriftPost, 'o', 'Color', 'b');
            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end), 'Pre', 'Color', 'r');
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end), 'Post', 'Color', 'b');
        else
            plot(ax_drift, crDriftPre, 'o', 'Color', MRSCont.colormap.Foreground);
            plot(ax_drift, crDriftPost, 'o', 'Color', MRSCont.colormap.Foreground, 'MarkerFaceColor', MRSCont.colormap.Foreground);
            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end), 'Pre', 'Color', MRSCont.colormap.Foreground);
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end), 'Post', 'Color', MRSCont.colormap.Foreground);
        end
        set(ax_drift, 'YLim', [3.028-0.06 3.028+0.06]);
        x = xlim;
        if ~GUI 
            plot(ax_drift, [x(1) x(2)], [3.028 3.028], ':k', 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04], '--k', 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04], '--k', 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        else
            plot(ax_drift, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color',MRSCont.colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color',MRSCont.colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color',MRSCont.colormap.Foreground, 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        end
    else
        axes(ax_drift);
        x = xlim;
        y = ylim;
        if ~GUI 
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color');            
        else
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color',MRSCont.colormap.Foreground);
        end
    end
else
    axes(ax_drift);
    x = xlim;
    y = ylim;
    if ~GUI 
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color');            
    else
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color',MRSCont.colormap.Foreground);
    end
end
    if ~GUI
        xlabel(ax_drift, 'Averages');
        ylabel(ax_drift, 'Cr frequency (ppm)');
        title(ax_drift, 'Frequency drift');
    else
        xlabel(ax_drift, 'Averages', 'Color', MRSCont.colormap.Foreground);
        ylabel(ax_drift, 'Cr frequency (ppm)', 'Color', MRSCont.colormap.Foreground);
        title(ax_drift, 'Frequency drift', 'Color', MRSCont.colormap.Foreground);        
    end


%%% 8. DESIGN FINETUNING %%%
% Adapt common style for all axes
axs = {ax_raw, ax_aligned, ax_proc, ax_drift};
for ll = 1:length(axs)
    gca = axs{ll};
    set(gca, 'LineWidth', 1, 'TickDir', 'out');
    set(gca, 'FontSize', 16);

    % Black axes, white background
    if ~GUI
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        % If no y caption, remove y axis
        if isempty(gca.YLabel.String)
            set(gca, 'YColor', 'w');
        else
            set(gca, 'YColor', 'k');
        end
    else
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        % If no y caption, remove y axis
        if isempty(gca.YLabel.String)
            set(gca, 'YColor', MRSCont.colormap.Background);
        else
            set(gca, 'YColor', MRSCont.colormap.Foreground);
        end        
    end

end

gcf = out;
    if ~GUI
        set(gcf, 'Color', 'w');
    else
        set(gcf, 'Color', MRSCont.colormap.Background);        
    end
box off;


%%% 9. ADD OSPREY LOGO %%%
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   