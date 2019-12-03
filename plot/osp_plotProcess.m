function out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax)
%% out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax)
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
if nargin<5
    switch which
        case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
            ppmmax = 4.5;
        case {'ref', 'w'}
            ppmmax = 2*4.68;
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    if nargin<4
        switch which
            case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
                ppmmin = 0.2;
            case {'ref', 'w'}
                ppmmin = 0;
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end
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


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range
switch which
    case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
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
out = figure;
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
for rr = 1:nAvgs
    plot(ax_raw, rawDataToPlot.ppm, rawDataToPlot.specs(:,rr), 'LineWidth', 0.5);
end
plotRange = op_freqrange(rawDataToPlot, ppmmin, ppmmax);
yLims = [mean(min(real(plotRange.specs))) mean(max(real(plotRange.specs)))];
set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
y = ylim;
plot(ax_raw, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_raw, [3.027 3.027], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_raw, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
hold(ax_raw, 'off');
title(ax_raw, 'Pre-alignment');


%%% 5. PLOT RAW ALIGNED %%%
% Apply stored corrections to calculate the spectra to display
applyDataToPlot = rawDataToPlot;
t = rawDataToPlot.t;
fs = procDataToPlot.specReg.fs;
phs = procDataToPlot.specReg.phs;
refShift = -repmat(MRSCont.QM.freqShift.(which)(kk), size(fs));
fs = fs + refShift;
for jj = 1:size(applyDataToPlot.fids,2)
    applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
end
applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);

hold(ax_aligned, 'on');    
% Loop over all averages
for rr = 1:nAvgs
    plot(ax_aligned, applyDataToPlot.ppm, applyDataToPlot.specs(:,rr), 'LineWidth', 0.5);
end
set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
y = ylim;
plot(ax_aligned, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_aligned, [3.027 3.027], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_aligned, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
hold(ax_aligned, 'off');
title(ax_aligned, 'Post-alignment');


%%% 6. PLOT PROCESSED %%%
% Add the data and plot
hold(ax_proc, 'on');    
plot(ax_proc, procDataToPlot.ppm, procDataToPlot.specs, 'k', 'LineWidth', 1.5);
set(ax_proc, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
y = ylim;
plot(ax_proc, [2.008 2.008], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_proc, [3.028 3.028], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
plot(ax_proc, [3.200 3.200], [y(1)-y(2) y(2)], ':k', 'LineWidth', 0.5);
hold(ax_proc, 'off');
title(ax_proc, 'Aligned and averaged');


%%% 7. GENERATE DRIFT PLOT %%%
if isfield(MRSCont.QM, 'drift')
    crDriftPre = MRSCont.QM.drift.(which)(:,kk);
    crDriftPost = op_measureDrift(applyDataToPlot);
    hold(ax_drift, 'on');
    plot(ax_drift, crDriftPre, 'o', 'Color', 'r');
    plot(ax_drift, crDriftPost, 'o', 'Color', 'b');
    text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end), 'Pre', 'Color', 'r');
    text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end), 'Post', 'Color', 'b');
    set(ax_drift, 'YLim', [3.028-0.06 3.028+0.06]);
    x = xlim;
    plot(ax_drift, [x(1) x(2)], [3.028 3.028], ':k', 'LineWidth', 0.5);
    plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04], '--k', 'LineWidth', 0.5);
    plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04], '--k', 'LineWidth', 0.5);
    hold(ax_drift, 'off');
else
    axes(ax_drift);
    x = xlim;
    y = ylim;
    text(ax_drift, x/2, y/2, 'No drift data available');
end
xlabel(ax_drift, 'Averages');
ylabel(ax_drift, 'Cr frequency (ppm)');
title(ax_drift, 'Frequency drift');


%%% 8. DESIGN FINETUNING %%%
% Adapt common style for all axes
axs = {ax_raw, ax_aligned, ax_proc, ax_drift};
for ll = 1:length(axs)
    gca = axs{ll};
    set(gca, 'LineWidth', 1, 'TickDir', 'out');
    set(gca, 'FontSize', 16);
    % If no y caption, remove y axis
    if isempty(gca.YLabel.String)
        set(gca, 'YColor', 'w');
    else
        set(gca, 'YColor', 'k');
    end
    % Black axes, white background
    set(gca, 'XColor', 'k');
    set(gca, 'Color', 'w');
    
end

gcf = out;
set(gcf, 'Color', 'w');
box off;


%%% 9. ADD OSPREY LOGO %%%
[I, map] = imread('osprey.gif','gif');
axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
imshow(I, map);
axis off;

end

   