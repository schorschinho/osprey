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
        case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
            ppmmax = 4.5;
        case {'ref', 'w'}
            ppmmax = 2*4.68;
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    if nargin<5
        switch which
            case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
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

% Set up colormaps
if GUI
    colormap = MRSCont.colormap;
else
    colormap.Background     = [1 1 1];
    colormap.LightAccent    = [110/255 136/255 164/255];
    colormap.Foreground     = [0 0 0];
    colormap.Accent         = [11/255 71/255 111/255];
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range
switch which
    case {'A', 'B', 'C', 'D'}
        raw            = MRSCont.raw{kk};
        procDataToPlot = MRSCont.processed.(which){kk};
        
        % Get sub-spectra, depending on whether they are stored as such
        if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
            if raw.subspecs == 4
                raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
                raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
                raw_C   = op_takesubspec(raw,3);                    % Get third subspectrum
                raw_D   = op_takesubspec(raw,4);                    % Get fourth subspectrum
            else
                raw_A   = op_takeaverages(raw,1:4:raw.averages);    % Get first subspectrum
                raw_B   = op_takeaverages(raw,2:4:raw.averages);    % Get second subspectrum
                raw_C   = op_takeaverages(raw,3:4:raw.averages);    % Get third subspectrum
                raw_D   = op_takeaverages(raw,4:4:raw.averages);    % Get fourth subspectrum
            end
        elseif MRSCont.flags.isMEGA
            if raw.subspecs == 2
                raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
                raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
            else
                raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
                raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
            end
        end
        
        eval(['rawDataToPlot = raw_' which ';']);
        
    case {'diff1', 'diff2', 'sum'}
        rawDataToPlot  = MRSCont.raw{kk};
        procDataToPlot = MRSCont.processed.(which){kk};
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

% Divide the figure into six tiles, create four axes
ax_raw      = subplot(3, 2, 1);
ax_aligned  = subplot(3, 2, 3);
ax_proc     = subplot(3, 2, 5);
ax_drift    = subplot(3, 2, 2);


%%% 4. PLOT RAW UNALIGNED %%%
% Add the data and plot
hold(ax_raw, 'on');    
% Loop over all averages
nAvgsRaw = rawDataToPlot.averages;
if ~GUI
    for rr = 1:nAvgsRaw
        plot(ax_raw, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5);
    end
else
    for rr = 1:nAvgsRaw
        plot(ax_raw, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5, 'Color', colormap.Foreground);
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
        plot(ax_raw, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
        plot(ax_raw, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
        plot(ax_raw, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
    end
        
end
hold(ax_raw, 'off');
title(ax_raw, 'Pre-alignment', 'Color', colormap.Foreground);
xlabel(ax_raw, 'Frequency (ppm)', 'Color', colormap.Foreground)
if GUI
    set(ax_raw, 'YColor', colormap.Background);
    set(ax_raw,'YTickLabel',{})
    set(ax_raw,'YTick',{})
end


%%% 5. PLOT RAW ALIGNED %%%
% Apply stored corrections to calculate the spectra to display
applyDataToPlot = rawDataToPlot;
t = rawDataToPlot.t;
switch which
    case {'A', 'B', 'C', 'D'} 
        fs = procDataToPlot.specReg.fs;
        phs = procDataToPlot.specReg.phs;
    case {'diff1', 'diff2', 'sum'}
        fs = rawDataToPlot.specReg.fs;
        phs = rawDataToPlot.specReg.phs;
end

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
    for rr = 1:nAvgsRaw
        plot(ax_aligned, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5);
    end
else
    for rr = 1:nAvgsRaw
        plot(ax_aligned, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5, 'Color', colormap.Foreground);
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
        plot(ax_aligned, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_aligned, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_aligned, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);        
    end
end
hold(ax_aligned, 'off');
title(ax_aligned, 'Post-alignment', 'Color', colormap.Foreground);
xlabel(ax_aligned, 'Frequency (ppm)', 'Color', colormap.Foreground)
if GUI
    set(ax_aligned, 'YColor', colormap.Background);
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
        plot(ax_proc, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_proc, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
        plot(ax_proc, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);        
    end    
end
hold(ax_proc, 'off');
title(ax_proc, 'Aligned and averaged', 'Color', colormap.Foreground);
xlabel(ax_proc, 'Frequency (ppm)', 'Color', colormap.Foreground)
if GUI
    set(ax_proc, 'YColor', colormap.Background);
    set(ax_proc,'YTickLabel',{})
    set(ax_proc,'YTick',{})
end

%%% 7. GENERATE DRIFT PLOT %%%
if isfield(MRSCont.QM.drift, which)
    if length(MRSCont.QM.drift.(which){kk}) > 1
        crDriftPre = MRSCont.QM.drift.pre.(which){kk};
        crDriftPost = MRSCont.QM.drift.post.(which){kk};
        hold(ax_drift, 'on');
        if ~GUI
            plot(ax_drift, crDriftPre, 'o', 'Color', 'r');
            plot(ax_drift, crDriftPost, 'o', 'Color', 'b');
            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end), 'Pre', 'Color', 'r');
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end), 'Post', 'Color', 'b');
        else
            plot(ax_drift, crDriftPre, 'o', 'Color', colormap.Foreground);
            plot(ax_drift, crDriftPost, 'o', 'Color', colormap.Foreground, 'MarkerFaceColor', colormap.Foreground);
            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end), 'Pre', 'Color', colormap.Foreground);
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end), 'Post', 'Color', colormap.Foreground);
        end
        set(ax_drift, 'YLim', [3.028-0.06 3.028+0.06]);
        x = xlim;
        if ~GUI 
            plot(ax_drift, [x(1) x(2)], [3.028 3.028], ':k', 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04], '--k', 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04], '--k', 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        else
            plot(ax_drift, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        end
    else
        axes(ax_drift);
        x = xlim;
        y = ylim;
        if ~GUI 
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color');            
        else
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
        end
    end
else
    axes(ax_drift);
    x = xlim;
    y = ylim;
    if ~GUI 
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color');            
    else
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
    end
end
    if ~GUI
        xlabel(ax_drift, 'Averages');
        ylabel(ax_drift, 'Cr frequency (ppm)');
        title(ax_drift, 'Frequency drift');
    else
        xlabel(ax_drift, 'Averages', 'Color', colormap.Foreground);
        ylabel(ax_drift, 'Cr frequency (ppm)', 'Color', colormap.Foreground);
        title(ax_drift, 'Frequency drift', 'Color', colormap.Foreground);        
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
        set(gca, 'XColor', colormap.Foreground);
        set(gca, 'Color', colormap.Background);
        % If no y caption, remove y axis
        if isempty(gca.YLabel.String)
            set(gca, 'YColor', colormap.Background);
        else
            set(gca, 'YColor', colormap.Foreground);
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
% Add to the printout, but not if displayed in the GUI.
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   