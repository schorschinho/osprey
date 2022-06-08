function [MRSCont] = OspreyHTMLReport(MRSCont,kk)
%% [MRSCont] = OspreyHTMLReport(MRSCont,kk)
%   This function creates a short HTML report of the processing and modeling
%   and should be called at the end of the analysis.
%
%   USAGE:
%       MRSCont = OspreyHTMLReport(MRSCont,kk);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = subject index
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2022-05-17)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2022-05-17: First version of the code.

colormaps = MRSCont.colormap;
ppmmin = 0.2;
ppmmax=4.2;

outputFolder    = fullfile(MRSCont.outputFolder,'Reports');
split_subject_path = strsplit(fileparts(MRSCont.files{kk}),filesep);
str_ind_sub = find(contains(split_subject_path,'sub'));
if ~isempty(str_ind_sub)
    sub_str = split_subject_path(str_ind_sub(1));
    sub_str = sub_str{1};
else
    sub_str = ['sub-' num2str(kk)];
end
outputFigures   = fullfile(MRSCont.outputFolder,'Reports','reportFigures',sub_str);
[foldername,filename,~]  = fileparts(MRSCont.files{kk});

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if ~exist(outputFigures,'dir')
    mkdir(outputFigures);
end

%% Process images export
if MRSCont.processed.metab{kk}.flags.isUnEdited
    ExtraIndex = 1;
    rawDataToPlot  = MRSCont.raw{ExtraIndex,kk};
    which_spec = 'metab';
    which_sub_spec ='A';
    procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
    SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
    procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
    rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
    applyDataToPlot = rawDataToPlot;
    t = rawDataToPlot.t;

    fs = procDataToPlot.specReg{1}.fs(:,SubSpectraIndex);
    phs = procDataToPlot.specReg{1}.phs(:,SubSpectraIndex);  
    weights = MRSCont.processed.(which_spec){kk}.specReg{1}.weights{find(strcmp(which_sub_spec,{'A', 'B', 'C', 'D'}))};      

    refShift = -repmat(MRSCont.QM.freqShift.(which_spec)(ExtraIndex,kk,SubSpectraIndex), size(fs));
    fs = fs - refShift;
    for jj = 1:size(applyDataToPlot.fids,2)
        applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
            exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
    end

    applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);
    plotRangeScale = op_freqrange(applyDataToPlot, ppmmin, ppmmax);

    yLims= [ min(min(real(plotRangeScale.specs(:,:)))) max(max(real(plotRangeScale.specs(:,:))))];
    yLimsAbs = (abs(yLims(1)) +  abs(yLims(2)));

    yLims = [yLims(1) - (yLimsAbs*0.1) yLims(2) + (yLimsAbs*0.1)];
    
    % Add the data and plot
    out = figure('Visible','off');
    hold(gca, 'on');    
    % Loop over all averages
    try
        nAvgsRaw = rawDataToPlot.sz(rawDataToPlot.dims.averages);
    catch % This is a wild guess in case no averages dimension is stored 
        nAvgsRaw = rawDataToPlot.sz(2);
    end


    for rr = 1:nAvgsRaw
        plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
        plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);           
    end
    set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);
    
     hold(gca, 'off');
    title(gca, 'Post-alignment', 'Color', colormaps.Foreground);
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'YColor', colormaps.Background);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


     box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Aligned_',which_spec]),'jpg');
    close(out);

    
    out = figure('Visible','off');
    hold(gca, 'on');  
    
    if isfield(MRSCont.QM.drift.pre, which_sub_spec)
        if length(MRSCont.QM.drift.pre.(which_sub_spec){kk}) > 1
            crDriftPre = MRSCont.QM.drift.pre.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
            crDriftPost = MRSCont.QM.drift.post.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
            colors = ones(length(crDriftPre),1).*colormaps.Foreground;
            for dots = 1 : length(crDriftPre)
                colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots,1));
                colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots,1));
                colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots,1));
            end
            scatter(gca, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormaps.LightAccent);
            scatter(gca, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormaps.Foreground);    
        
            text(gca, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormaps.LightAccent);
            text(gca, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormaps.Foreground);
            set(gca, 'YLim', [3.028-0.1 3.028+0.1]);
            yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
            yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
            x = xlim;
            plot(gca, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
            plot(gca, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
            plot(gca, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
            hold(gca, 'off');
        else 
            x = xlim;
            y = yLims;
            text(gca, x(2)/6, y(2)/2, 'No drift data available','Color', colormaps.Foreground);
        end
    else
        x = xlim;
        y = yLims;
        text(gca, x(2)/6, y(2)/2, 'No drift data available','Color', colormaps.Foreground);
    end
    xlabel(gca, 'Averages', 'Color', colormaps.Foreground);
    ylabel(gca, 'Cr frequency (ppm)', 'Color', colormaps.Foreground);
    title(gca, 'chemical shift drift', 'Color', colormaps.Foreground); 

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'YColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


    box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Drift_',which_spec]),'jpg');
    close(out);

    ProcSpecNames = {'A'};
    for ss = 1 : length(ProcSpecNames)
        which_sub_spec = ProcSpecNames{ss};
        procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
        SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
        procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));        
        out = figure('Visible','off');
        hold(gca, 'on');    
        plot(gca, procDataToPlot.ppm, real(procDataToPlot.specs(:,1))/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
        y = [-0.2, 1.2];
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', y);
        if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'MM_ref'))
            plot(gca, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            plot(gca, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            if ~strcmp(which_sub_spec, 'mm')
                plot(gca, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5); 
            else
                plot(gca, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            end
        end
        hold(gca, 'off');
        title(gca, which_sub_spec, 'Color', colormaps.Foreground);
        xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
        set(gca, 'YColor', colormaps.Background);
        set(gca,'YTickLabel',{})
        set(gca,'YTick',{})
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
        set(gca, 'FontSize', 16);
        set(gca, 'XColor', colormaps.Foreground);
        set(gca, 'Color', colormaps.Background);
    
        box off;
        set(out,'PaperUnits','centimeters');
        set(out,'PaperPosition',[0 0 20 15]);
        saveas(out,fullfile(outputFigures,[sub_str '_' which_spec,'_',which_sub_spec]),'jpg');
        close(out);
    end
end

if MRSCont.processed.metab{kk}.flags.isMEGA
    ExtraIndex = 1;
    rawDataToPlot  = MRSCont.raw{ExtraIndex,kk};
    which_spec = 'metab';
    which_sub_spec ='sum';
    procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
    SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
    procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
    if procDataToPlot.flags.orderswitched
        temp_spec = rawDataToPlot.specs(:,:,1);
        rawDataToPlot.specs(:,:,1) = rawDataToPlot.specs(:,:,2);
        rawDataToPlot.specs(:,:,2) = temp_spec;
        temp_fids = rawDataToPlot.fids(:,:,1);
        rawDataToPlot.fids(:,:,1) = rawDataToPlot.fids(:,:,2);
        rawDataToPlot.fids(:,:,2) = temp_fids;
    end
    proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                   % Get first subspectrum
    proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                  % Get second subspectrum
    rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
    applyDataToPlot = rawDataToPlot;
    t = rawDataToPlot.t;

    fs{1} = proc_A.specReg{1}.fs;
    phs{1} = proc_A.specReg{1}.phs;
    fs{2} = proc_B.specReg{1}.fs;
    phs{2} = proc_B.specReg{1}.phs;
    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{1}.weights{:});

    refShift = -repmat(MRSCont.QM.freqShift.(which_spec)(ExtraIndex,kk,SubSpectraIndex), size(fs{1}));
    for ss = 1 : length(fs)
        fs{ss} = fs{ss} - refShift;
        for jj = 1:size(applyDataToPlot.fids,2)
            if length(size(applyDataToPlot.fids)) == 3
                applyDataToPlot.fids(:,jj,ss) = applyDataToPlot.fids(:,jj,ss) .* ...
                    exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
            else
                applyDataToPlot.fids(:,ss) = applyDataToPlot.fids(:,ss) .* ...
                    exp(1i*fs{ss}(1)*2*pi*t') * exp(1i*pi/180*phs{ss}(1));
            end
        end
    end
    applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);

    plotRangeScale = op_freqrange(applyDataToPlot, ppmmin, ppmmax);

    yLims= [ min(min(real(plotRangeScale.specs(:,:)))) max(max(real(plotRangeScale.specs(:,:))))];
    yLimsAbs = (abs(yLims(1)) +  abs(yLims(2)));
    if strcmp(which_sub_spec, 'diff1')
        yLims = [yLims(1) - (yLimsAbs*0.1) (2*yLims(2)) + (yLimsAbs*0.1)];
    else
        yLims = [yLims(1) - (yLimsAbs*0.1) yLims(2) + (yLimsAbs*0.1)];
    end
    
    % Add the data and plot
    out = figure('Visible','off');
    hold(gca, 'on');    
    % Loop over all averages
    try
        nAvgsRaw = rawDataToPlot.sz(rawDataToPlot.dims.averages);
    catch % This is a wild guess in case no averages dimension is stored 
        nAvgsRaw = rawDataToPlot.sz(2);
    end


     if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref') && ~strcmp(which_sub_spec, 'A')
        stag = [0,0.5,1,1.5] .* yLimsAbs;
        stagText = stag + (0.25.* yLimsAbs);
        for rr = 1:nAvgsRaw
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormaps.LightAccent);
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
        end
        text(gca, ppmmin+0.3, stagText(1), 'A', 'Color', colormaps.LightAccent);
        text(gca, ppmmin+0.3, stagText(2), 'B', 'Color', colormaps.Foreground);
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);           
    else
        for rr = 1:nAvgsRaw
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
            if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref')
                plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);           
            end
        end
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);
     end
    
     hold(gca, 'off');
    title(gca, 'Post-alignment', 'Color', colormaps.Foreground);
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'YColor', colormaps.Background);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


     box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Aligned_',which_spec]),'jpg');
    close(out);

    
    out = figure('Visible','off');
    hold(gca, 'on');  

    crDriftPre = MRSCont.QM.drift.pre.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
    crDriftPost = MRSCont.QM.drift.post.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
    colors = ones(length(crDriftPre),1).*colormaps.Foreground;
    for dots = 1 : length(crDriftPre)
        colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots,1));
        colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots,1));
        colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots,1));
    end
    scatter(gca, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormaps.LightAccent);
    scatter(gca, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormaps.Foreground);    

    text(gca, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormaps.LightAccent);
    text(gca, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormaps.Foreground);
    set(gca, 'YLim', [3.028-0.1 3.028+0.1]);
    yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
    yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
    x = xlim;
    plot(gca, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    plot(gca, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    plot(gca, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    hold(gca, 'off');

    xlabel(gca, 'Averages', 'Color', colormaps.Foreground);
    ylabel(gca, 'Cr frequency (ppm)', 'Color', colormaps.Foreground);
    title(gca, 'chemical shift drift', 'Color', colormaps.Foreground); 

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'YColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


    box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Drift_',which_spec]),'jpg');
    close(out);

    ProcSpecNames = {'A','diff1'};
    for ss = 1 : length(ProcSpecNames)
        which_sub_spec = ProcSpecNames{ss};
        procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
        SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
        procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));        
        out = figure('Visible','off');
        hold(gca, 'on');    
        plot(gca, procDataToPlot.ppm, real(procDataToPlot.specs(:,1))/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
        if strcmp(which_sub_spec,'diff2')
            y = [-1.2, 1.2];
        else if strcmp(which_sub_spec,'diff1')
                if strcmp(which_spec,'metab')
                    y = [-2, 1.2];   
                else
                    y = [-1.5, 1.2];   
                end
            else
                if strcmp(which_spec,'metab')
                    y = [-0.2, 1.2];
                else
                    y = [-1.5, 1.2];   
                end
    
            end
        end
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', y);
        if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'MM_ref'))
            plot(gca, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            plot(gca, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            if ~strcmp(which_sub_spec, 'mm')
                plot(gca, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5); 
            else
                plot(gca, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            end
        end
        hold(gca, 'off');
        title(gca, which_sub_spec, 'Color', colormaps.Foreground);
        xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
        set(gca, 'YColor', colormaps.Background);
        set(gca,'YTickLabel',{})
        set(gca,'YTick',{})
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
        set(gca, 'FontSize', 16);
        set(gca, 'XColor', colormaps.Foreground);
        set(gca, 'Color', colormaps.Background);
    
        box off;
        set(out,'PaperUnits','centimeters');
        set(out,'PaperPosition',[0 0 20 15]);
        saveas(out,fullfile(outputFigures,[sub_str '_' which_spec '_' which_sub_spec]),'jpg');
        close(out);
    end
end

if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    ExtraIndex = 1;
    rawDataToPlot  = MRSCont.raw{ExtraIndex,kk};
    which_spec = 'metab';
    which_sub_spec ='sum';
    procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
    SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
    procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
    temp_spec = cat(3,rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(2)),...
                    rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(4)));
    temp_fid = cat(3,rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(2)),...
                    rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(4)));
    rawDataToPlot.fids = temp_fid;
    rawDataToPlot.specs = temp_spec;
    proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                   % Get first subspectrum
    proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                  % Get second subspectrum
    proc_C   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'C')));                   % Get third subspectrum
    proc_D   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'D')));                     % Get fourth subspectrum
    rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
    applyDataToPlot = rawDataToPlot;
    t = rawDataToPlot.t;

    fs{1} = proc_A.specReg{1}.fs;
    phs{1} = proc_A.specReg{1}.phs;
    fs{2} = proc_B.specReg{1}.fs;
    phs{2} = proc_B.specReg{1}.phs;
    fs{3} = proc_C.specReg{1}.fs;
    phs{3} = proc_C.specReg{1}.phs;
    fs{4} = proc_D.specReg{1}.fs;
    phs{4} = proc_D.specReg{1}.phs;
    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{1}.weights{:});

    refShift = -repmat(MRSCont.QM.freqShift.(which_spec)(ExtraIndex,kk,SubSpectraIndex), size(fs{1}));
    for ss = 1 : length(fs)
        fs{ss} = fs{ss} - refShift;
        for jj = 1:size(applyDataToPlot.fids,2)
            if length(size(applyDataToPlot.fids)) == 3
                applyDataToPlot.fids(:,jj,ss) = applyDataToPlot.fids(:,jj,ss) .* ...
                    exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
            else
                applyDataToPlot.fids(:,ss) = applyDataToPlot.fids(:,ss) .* ...
                    exp(1i*fs{ss}(1)*2*pi*t') * exp(1i*pi/180*phs{ss}(1));
            end
        end
    end
    applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);

    plotRangeScale = op_freqrange(applyDataToPlot, ppmmin, ppmmax);

    yLims= [ min(min(real(plotRangeScale.specs(:,:)))) max(max(real(plotRangeScale.specs(:,:))))];
    yLimsAbs = (abs(yLims(1)) +  abs(yLims(2)));
    if strcmp(which_sub_spec, 'diff1') || strcmp(which_sub_spec, 'diff2') || strcmp(which_sub_spec, 'diff3') || strcmp(which_sub_spec, 'sum')
        yLims = [yLims(1) - (yLimsAbs*0.1) (3*yLims(2)) + (yLimsAbs*0.1)];
    else
        yLims = [yLims(1) - (yLimsAbs*0.1) yLims(2) + (yLimsAbs*0.1)];
    end
    
    % Add the data and plot
    out = figure('Visible','off');
    hold(gca, 'on');    
    % Loop over all averages
    try
        nAvgsRaw = rawDataToPlot.sz(rawDataToPlot.dims.averages);
    catch % This is a wild guess in case no averages dimension is stored 
        nAvgsRaw = rawDataToPlot.sz(2);
    end


     if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref') && ~strcmp(which_sub_spec, 'A') && ~strcmp(which_sub_spec, 'B') && ~strcmp(which_sub_spec, 'C') && ~strcmp(which_sub_spec, 'D')
        stag = [0,0.5,1,1.5] .* yLimsAbs;
        stagText = stag + (0.25.* yLimsAbs);
        for rr = 1:nAvgsRaw
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormaps.LightAccent);
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,3) + stag(3)), 'LineWidth', 0.5, 'Color', colormaps.LightAccent);
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,4) + stag(4)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
        end
        text(gca, ppmmin+0.3, stagText(1), 'A', 'Color', colormaps.LightAccent);
        text(gca, ppmmin+0.3, stagText(2), 'B', 'Color', colormaps.Foreground);
        text(gca, ppmmin+0.3, stagText(3), 'C', 'Color', colormaps.LightAccent);
        text(gca, ppmmin+0.3, stagText(4), 'D', 'Color', colormaps.Foreground); 
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);           
    else
        for rr = 1:nAvgsRaw
            plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);
            if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref')
                plot(gca, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormaps.Foreground);           
            end
        end
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);
     end
    
     hold(gca, 'off');
    title(gca, 'Post-alignment', 'Color', colormaps.Foreground);
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'YColor', colormaps.Background);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


     box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Aligned_',which_spec]),'jpg');
    close(out);

    
    out = figure('Visible','off');
    hold(gca, 'on');  

    crDriftPre = MRSCont.QM.drift.pre.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
    crDriftPost = MRSCont.QM.drift.post.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
    colors = ones(length(crDriftPre),1).*colormaps.Foreground;
    for dots = 1 : length(crDriftPre)
        colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots,1));
        colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots,1));
        colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots,1));
    end
    scatter(gca, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormaps.LightAccent);
    scatter(gca, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormaps.Foreground);    

    text(gca, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormaps.LightAccent);
    text(gca, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormaps.Foreground);
    set(gca, 'YLim', [3.028-0.1 3.028+0.1]);
    yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
    yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
    x = xlim;
    plot(gca, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    plot(gca, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    plot(gca, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormaps.Foreground, 'LineWidth', 0.5);
    hold(gca, 'off');

    xlabel(gca, 'Averages', 'Color', colormaps.Foreground);
    ylabel(gca, 'Cr frequency (ppm)', 'Color', colormaps.Foreground);
    title(gca, 'chemical shift drift', 'Color', colormaps.Foreground); 

    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'YColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);


    box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_Drift_',which_spec]),'jpg');
    close(out);

    ProcSpecNames = {'sum','diff1','diff2'};
    for ss = 1 : length(ProcSpecNames)
        which_sub_spec = ProcSpecNames{ss};
        procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
        SubSpectraIndex = find(strcmp(which_sub_spec,procDataToPlot.names));
        procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));        
        out = figure('Visible','off');
        hold(gca, 'on');    
        plot(gca, procDataToPlot.ppm, real(procDataToPlot.specs(:,1))/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
        if strcmp(which_sub_spec,'diff2')
            y = [-1.2, 1.2];
        else if strcmp(which_sub_spec,'diff1')
                if strcmp(which_spec,'metab')
                    y = [-2, 1.2];   
                else
                    y = [-1.5, 1.2];   
                end
            else
                if strcmp(which_spec,'metab')
                    y = [-0.2, 1.2];
                else
                    y = [-1.5, 1.2];   
                end
    
            end
        end
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', y);
        if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'MM_ref'))
            plot(gca, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            plot(gca, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            if ~strcmp(which_sub_spec, 'mm')
                plot(gca, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5); 
            else
                plot(gca, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormaps.Foreground,  'LineWidth', 0.5);
            end
        end
        hold(gca, 'off');
        title(gca, which_sub_spec, 'Color', colormaps.Foreground);
        xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
        set(gca, 'YColor', colormaps.Background);
        set(gca,'YTickLabel',{})
        set(gca,'YTick',{})
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
        set(gca, 'FontSize', 16);
        set(gca, 'XColor', colormaps.Foreground);
        set(gca, 'Color', colormaps.Background);
    
        box off;
        set(out,'PaperUnits','centimeters');
        set(out,'PaperPosition',[0 0 20 15]);
        saveas(out,fullfile(outputFigures,[sub_str '_' which_spec,'_',which_sub_spec]),'jpg');
        close(out);
    end
end

if MRSCont.flags.hasRef
    ppmmin = 0;
    ppmmax = 2*4.68;
    which_spec = 'ref';
    procDataToPlot = op_takeextra(MRSCont.processed.ref{kk},ExtraIndex); 
    out = figure('Visible','off');
    hold(gca, 'on');   

    plot(gca, procDataToPlot.ppm, real(procDataToPlot.specs(:,1)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
    y = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)];
    set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);

    hold(gca, 'off');
    title(gca, 'eddy current reference', 'Color', colormaps.Foreground);
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'YColor', colormaps.Background);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})
    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);

    box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_' which_spec]),'jpg');
    close(out);

end

if MRSCont.flags.hasWater
    pmmin = 0;
    ppmmax = 2*4.68;
    which_spec = 'w';
    procDataToPlot = op_takeextra(MRSCont.processed.ref{kk},ExtraIndex); 
    out = figure('Visible','off');
    hold(gca, 'on');   

    plot(gca, procDataToPlot.ppm, real(procDataToPlot.specs(:,1)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
    y = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)];
    set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);

    hold(gca, 'off');
    title(gca, 'short TE water', 'Color', colormaps.Foreground);
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'YColor', colormaps.Background);
    set(gca,'YTickLabel',{})
    set(gca,'YTick',{})
    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);
    set(gca, 'XColor', colormaps.Foreground);
    set(gca, 'Color', colormaps.Background);

    box off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_' which_spec]),'jpg');
    close(out);

end
%% Fit images export
% model_names = {'sum','diff1','diff2'};
if MRSCont.processed.metab{kk}.flags.isUnEdited
    spec_names = {'metab'};
    VoxelIndices = {[1,1,1]};
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    spec_names = {'metab','metab'};
    VoxelIndices = {[1,1,1],...
                    [1,2,1]};
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    spec_names = {'metab','metab','metab'};
    VoxelIndices = {[1,1,1],...
                    [1,2,1],...
                    [1,3,1]};
end
if MRSCont.flags.hasRef
    spec_names{end+1}= 'ref';
    VoxelIndices{end+1}=[1,1,1];
end
if MRSCont.flags.hasWater
    spec_names{end+1}= 'w';
    VoxelIndices{end+1}=[1,1,1];
end
for f = 1 : length(spec_names)
    which_spec=spec_names{f};
    VoxelIndex = VoxelIndices{f};
    dataToPlot = MRSCont.processed.(which_spec){kk};
     if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
        fitRangePPM = MRSCont.opts.fit.rangeWater;
        basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]){VoxelIndex,1};
        fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),kk};
     else
        basisSet    = MRSCont.fit.resBasisSet.(which_spec).(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]){1,1,VoxelIndex(2)};
        dataToPlot   = op_takesubspec(MRSCont.processed.(which_spec){kk},find(strcmp(MRSCont.processed.(which_spec){kk}.names,basisSet.names{1})));
        fitParams   = MRSCont.fit.results.(which_spec).fitParams{VoxelIndex(3),kk,VoxelIndex(2)};
        fitRangePPM = MRSCont.opts.fit.range;
     end
    
    
    
    inputData.dataToFit                 = dataToPlot;
    inputData.basisSet                  = basisSet;
    inputSettings.scale                 = MRSCont.fit.scale{kk};
    inputSettings.fitRangePPM           = fitRangePPM;
    inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
    inputSettings.fitStyle              = MRSCont.opts.fit.style;
    inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
    inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
    inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
    inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
    inputSettings.concatenated.Subspec  = basisSet.names{1};
    if isfield(MRSCont.opts.fit,'GAP') && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
        inputSettings.GAP = MRSCont.opts.fit.GAP.(dataToPlot.names{1});
    else
        inputSettings.GAP = [];
    end
    if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
        [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
    else
        [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
    end
    nBasisFct = basisSet.nMets;
    if isfield(basisSet, 'nMM')
        nBasisFct = nBasisFct + basisSet.nMM;
    end
    colorData = MRSCont.colormap.Foreground;
    colorFit  = MRSCont.colormap.Accent;
    linewidthFit = 1.6;
    linewidthResidual = 1;
    colorBaseline = MRSCont.colormap.LightAccent;
    ppm         = ModelOutput.ppm;
    dataToPlot  = ModelOutput.data;
    fit         = ModelOutput.completeFit;
    residual    = ModelOutput.residual;
    if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
        baseline    = ModelOutput.baseline;
        indivPlots  = ModelOutput.indivMets;
        basisSetNames = basisSet.name;
    end
    stagData = 0.1*(max(abs(min(dataToPlot)), abs(max(dataToPlot))));
    maxPlot = max(dataToPlot + abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData;
    out = figure('Visible','off');
    hold(gca, 'on');  
    plot(ppm, (zeros(1,length(ppm)) + stagData)/maxPlot, 'Color', colorData); % Zeroline
    plot(ppm, (dataToPlot + stagData)/maxPlot, 'Color', colorData); % Data
    plot(ppm, (fit + stagData)/maxPlot, 'Color', colorFit, 'LineWidth', linewidthFit); % Fit
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', 1); % Maximum Data
    
    plot(ppm, (residual + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', linewidthResidual); % Residual
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color', colorData, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot, 'Color', colorData, 'LineWidth', 1); % Max Residue
    
    
    
    text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Zeroline text
    text(fitRangePPM(1), (0 + max(dataToPlot) + stagData)/maxPlot-0.05, num2str(max(dataToPlot),'%1.2e'), 'FontSize', 10,'Color', MRSCont.colormap.Foreground); % Maximum Data Text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '0', 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Zeroline Residual text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot +0.05, num2str(abs(max(dataToPlot - fit)),'%1.2e'), 'FontSize', 10,'Color', MRSCont.colormap.Foreground); %Max Residue text
    
    if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w') || contains(which_spec, 'mm')) 
        plot(ppm, (real(baseline) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', colorBaseline);   
        stag = maxPlot *  2.5 / nBasisFct;
        % Loop over all basis functions
        for rr = 1:nBasisFct
            plot(ppm, (indivPlots(:,rr) - rr*stag)/maxPlot, 'Color', MRSCont.colormap.Foreground);
            text(fitRangePPM(1), (- rr*stag)/maxPlot, basisSetNames{rr}, 'FontSize', 10,'Color', MRSCont.colormap.Foreground);
        end
    end
    % Adapt common style for all axes
    set(gca, 'XDir', 'reverse', 'XLim', [fitRangePPM(1), fitRangePPM(end)], 'XMinorTick', 'On');
    if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
        set(gca, 'YLim', [(-nBasisFct-1)*stag/maxPlot  1]);
    else
        set(gca, 'YLim', [0  1.2]);
    end
    
    set(gca, 'LineWidth', 1, 'TickDir', 'out');
    set(gca, 'FontSize', 16);
    set(gca, 'YColor', MRSCont.colormap.Background);
    set(gca,'YTickLabel',{});
    set(gca,'YTick',{});
    % Dirtywhite axes, light gray background
    xlabel(gca, 'chemical shift (ppm)', 'Color', colormaps.Foreground)
    set(gca, 'XColor', MRSCont.colormap.Foreground);
    set(gca, 'Color', MRSCont.colormap.Background);
    set(gcf, 'Color', MRSCont.colormap.Background);
    title(['Model ' which_spec ' ' basisSet.names{1}], 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
    box off;
    hold off;
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 20 15]);
    saveas(out,fullfile(outputFigures,[sub_str '_' which_spec '_' basisSet.names{1} '_model']),'jpg');
    close(out);
end
%% Coreg/Seg images export
if MRSCont.flags.didCoreg
    out = figure('Visible','off');
    MRSCont.flags.isGUI = 1;
    temp = osp_plotCoreg(MRSCont, kk);
    ViewAxes = gca();
    ViewAxes.Title.String ='';
    set(ViewAxes, 'Parent', out );
    colormap(out.Children,'gray');
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 18 6]);
    saveas(out,fullfile(outputFigures,[sub_str '_coreg_svs_space-scanner_mask']),'jpg');
    close(out);
    close(temp);
end
if MRSCont.flags.didSeg
    out = figure('Visible','off');
    temp = osp_plotSegment(MRSCont, kk);
    ViewAxes = gca();
    ViewAxes.Title.String ='';
    set(ViewAxes, 'Parent', out );
    colormap(out.Children,'gray');
    set(out,'PaperUnits','centimeters');
    set(out,'PaperPosition',[0 0 18 6]);
    saveas(out,fullfile(outputFigures,[sub_str '_seg_svs_space-scanner_mask']),'jpg');
    close(out);
    close(temp);
end
%% Write report in HTML files
%write an html report: 
fid=fopen(fullfile(outputFolder,[sub_str,'-report.html']),'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
fprintf(fid,'\n<head>');
fprintf(fid,'\n<style>');
fprintf(fid,'\n* {');
fprintf(fid,'\n\tbox-sizing: border-box;');
fprintf(fid,'\n}');
fprintf(fid,'\n\n.column {');
fprintf(fid,'\n\tfloat: left;');
fprintf(fid,'\nwidth: 50%%;');
fprintf(fid,'\npadding: 5px;');
fprintf(fid,'\n}');
fprintf(fid,'\n\n.row::after {');
fprintf(fid,'\n\tcontent: "";');
fprintf(fid,'\n\tclear: both;');
fprintf(fid,'\n\tdisplay: table;');
fprintf(fid,'\n}');
fprintf(fid,'\n\n.column3 {');
fprintf(fid,'\n\tfloat: left;');
fprintf(fid,'\n\twidth: 33.33%%;');
fprintf(fid,'\n\tpadding: 5px;');
fprintf(fid,'\n}');

fprintf(fid,'\n\n.row::after {');
fprintf(fid,'\n\tcontent: "";');
fprintf(fid,'\n\tclear: both;');
fprintf(fid,'\n\tdisplay: table;');
fprintf(fid,'\n}');
fprintf(fid,'\n</style>');
fprintf(fid,'\n</head>');
fprintf(fid,'\n<body>');

logoPath=which('osprey.png');
if ~isempty(logoPath)
    fprintf(fid,'\n<img src= " %s " width="35" height="28"> <b> \tOsprey Analysis Report</b> ',logoPath);
else
    fprintf(fid,'\n<b> Osprey Analysis Report</b>');
end
fprintf(fid,'\n<p><b>DATE:</b> %s \t <b>FILENAME:</b> %s </p>',date,filename);
fprintf(fid,'\n<h2>Summary:</h2>');
fprintf(fid,'\n<div class="row">');
fprintf(fid,'\n\t<div class="column3">');

fprintf(fid,'\n<p><b>signal-to-noise tCr</b> \t%5.2f',table2array(MRSCont.QM.tables(kk,1))); 
if table2array(MRSCont.QM.tables(kk,2)) / MRSCont.processed.metab{kk}.txfrq*1e6 < 0.1 
    fprintf(fid,'\n<p  style="color:green;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end
if table2array(MRSCont.QM.tables(kk,2))/ MRSCont.processed.metab{kk}.txfrq*1e6 > 0.1 && table2array(MRSCont.QM.tables(kk,2)) / MRSCont.processed.metab{kk}.txfrq*1e6 < 0.15
    fprintf(fid,'\n<p style="color:orange;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end
if table2array(MRSCont.QM.tables(kk,2))/ MRSCont.processed.metab{kk}.txfrq*1e6 > 0.15
    fprintf(fid,'\n<p style="color:red;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end

if MRSCont.flags.hasRef
    if table2array(MRSCont.QM.tables(kk,3)) / MRSCont.processed.ref{kk}.txfrq*1e6 < 0.1 
        fprintf(fid,'\n<p style="color:green;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3)));
    end
    if table2array(MRSCont.QM.tables(kk,3))/ MRSCont.processed.ref{kk}.txfrq*1e6 > 0.1 && table2array(MRSCont.QM.tables(kk,3)) / MRSCont.processed.ref{kk}.txfrq*1e6 < 0.15
        fprintf(fid,'\n<p style="color:orange;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3))/ MRSCont.processed.ref{kk}.txfrq*1e6);
    end   
    if table2array(MRSCont.QM.tables(kk,3))/ MRSCont.processed.ref{kk}.txfrq*1e6 > 0.15
        fprintf(fid,'\n<p style="color:red;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3)));
    end
end
fprintf(fid,'\n\t</div>'); 
fprintf(fid,'\n\t<div class="column3">');
if MRSCont.processed.metab{kk}.flags.isUnEdited
    fprintf(fid,'\n<p><b>Model Residual A [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    fprintf(fid,'\n<p><b>Model Residual diff1 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,7)));
    fprintf(fid,'\n<p><b>Model Residual A [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    fprintf(fid,'\n<p><b>Model Residual diff1 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,7)));
    fprintf(fid,'\n<p><b>Model Residual diff2 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,8)));
    fprintf(fid,'\n<p><b>Model Residual sum [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
fprintf(fid,'\n\t</div>');
fprintf(fid,'\n</div>');
fprintf(fid,'\n<div class="row">');
if MRSCont.processed.metab{kk}.flags.isUnEdited
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A_model.jpg']));
    fprintf(fid,'\n\t</div>');    
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A_model.jpg']));
    fprintf(fid,'\n\t</div>');
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff2_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_sum_model.jpg']));
    fprintf(fid,'\n\t</div>');
end
fprintf(fid,'\n</div>');
if MRSCont.flags.hasRef || MRSCont.flags.hasWater
    fprintf(fid,'\n<div class="row">');
    if MRSCont.flags.hasRef 
        fprintf(fid,'\n\t<div class="column3">');
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_ref_A_model.jpg']));
        fprintf(fid,'\n\t</div>');
        if MRSCont.flags.hasWater
            fprintf(fid,'\n\t<div class="column3">');
            fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_w_A_model.jpg']));
            fprintf(fid,'\n\t</div>');
        end
    
    end
if MRSCont.flags.didCoreg
    fprintf(fid,'\n\t<div class="column3">');
    if MRSCont.flags.didSeg
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%"> \n\t\t <img src= " %s" style="width:100%%">',...
            fullfile(outputFigures,[sub_str '_coreg_svs_space-scanner_mask.jpg']),...
            fullfile(outputFigures,[sub_str '_seg_svs_space-scanner_mask.jpg']));
    else
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_coreg_svs_space-scanner_mask.jpg']));
    end
    fprintf(fid,'\n\t</div>');
end
fprintf(fid,'\n</div>');
end
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<div class="row">');
fprintf(fid,'\n\t<div class="column3">');
fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_Aligned_metab.jpg']));
fprintf(fid,'\n\t</div>');
fprintf(fid,'\n\t<div class="column3">');
fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_Drift_metab.jpg']));
fprintf(fid,'\n\t</div>');
fprintf(fid,'\n</div>');
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Averaged spectra:</h2>');
fprintf(fid,'\n<p><b>signal-to-noise tCr</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,1))); 
if table2array(MRSCont.QM.tables(kk,2)) / MRSCont.processed.metab{kk}.txfrq*1e6 < 0.1 
    fprintf(fid,'\n<p style="color:green;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end
if table2array(MRSCont.QM.tables(kk,2))/ MRSCont.processed.metab{kk}.txfrq*1e6 > 0.1 && table2array(MRSCont.QM.tables(kk,2)) / MRSCont.processed.metab{kk}.txfrq*1e6 < 0.15
    fprintf(fid,'\n<p style="color:orange;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end
if table2array(MRSCont.QM.tables(kk,2))/ MRSCont.processed.metab{kk}.txfrq*1e6 > 0.15
    fprintf(fid,'\n<p style="color:red;"><b>linewidth tCr [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,2)));
end

if MRSCont.flags.hasRef
    if table2array(MRSCont.QM.tables(kk,3)) / MRSCont.processed.ref{kk}.txfrq*1e6 < 0.1 
        fprintf(fid,'\n<p style="color:green;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3)));
    end
    if table2array(MRSCont.QM.tables(kk,3))/ MRSCont.processed.ref{kk}.txfrq*1e6 > 0.1 && table2array(MRSCont.QM.tables(kk,3)) / MRSCont.processed.ref{kk}.txfrq*1e6 < 0.15
        fprintf(fid,'\n<p style="color:orange;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3)));
    end   
    if table2array(MRSCont.QM.tables(kk,3))/ MRSCont.processed.ref{kk}.txfrq*1e6 > 0.15
        fprintf(fid,'\n<p style="color:red;"><b>linewidth water [Hz]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,3)));
    end
end
if MRSCont.processed.metab{kk}.flags.isUnEdited
    fprintf(fid,'\n<div class="row">');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n</div>');
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    fprintf(fid,'\n<div class="row">');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n</div>');
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    fprintf(fid,'\n<div class="row">');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff2.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_sum.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n</div>');
end
if MRSCont.flags.hasRef || MRSCont.flags.hasWater
    fprintf(fid,'\n<div class="row">');
    if MRSCont.flags.hasRef
        fprintf(fid,'\n\t<div class="column3">');
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_ref.jpg']));
        fprintf(fid,'\n\t</div>');
    end
    if MRSCont.flags.hasWater
        fprintf(fid,'\n\t<div class="column3">');
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_w.jpg']));
        fprintf(fid,'\n\t</div>');
    end
    fprintf(fid,'\n</div>');
end
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Model results:</h2>');
if MRSCont.processed.metab{kk}.flags.isUnEdited
    fprintf(fid,'\n<p><b>Relative Model Residual A [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    fprintf(fid,'\n<p><b>Relative Model Residual diff1 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,7)));
    fprintf(fid,'\n<p><b>Relative Model Residual A [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    fprintf(fid,'\n<p><b>Relative Model Residual diff1 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,7)));
    fprintf(fid,'\n<p><b>Relative Model Residual diff2 [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,8)));
    fprintf(fid,'\n<p><b>Relative Model Residual sum [%%]</b> \t%5.2f </p>',table2array(MRSCont.QM.tables(kk,6)));
end
fprintf(fid,'\n<div class="row">');
if MRSCont.processed.metab{kk}.flags.isUnEdited
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A_model.jpg']));
    fprintf(fid,'\n\t</div>');    
end
if MRSCont.processed.metab{kk}.flags.isMEGA
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_A_model.jpg']));
    fprintf(fid,'\n\t</div>');
end
if MRSCont.processed.metab{kk}.flags.isHERMES || MRSCont.processed.metab{kk}.flags.isHERCULES
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff1_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_diff2_model.jpg']));
    fprintf(fid,'\n\t</div>');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_metab_sum_model.jpg']));
    fprintf(fid,'\n\t</div>');
end
fprintf(fid,'\n</div>');
if MRSCont.flags.hasRef || MRSCont.flags.hasWater
    fprintf(fid,'\n<div class="row">');
    if MRSCont.flags.hasRef 
        fprintf(fid,'\n\t<div class="column3">');
        fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_ref_A_model.jpg']));
        fprintf(fid,'\n\t</div>');
        if MRSCont.flags.hasWater
            fprintf(fid,'\n\t<div class="column3">');
            fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_w_A_model.jpg']));
            fprintf(fid,'\n\t</div>');
        end
    
    end
fprintf(fid,'\n</div>');
end
if MRSCont.flags.didCoreg
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<h2>Coregistration & Segmentation:</h2>');
    fprintf(fid,'\n<div class="row">');
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_coreg_svs_space-scanner_mask.jpg']));
    fprintf(fid,'\n\t</div>');
end
if MRSCont.flags.didSeg
    fprintf(fid,'\n\t<div class="column3">');
    fprintf(fid,'\n\t\t<img src= " %s" style="width:100%%">',fullfile(outputFigures,[sub_str '_seg_svs_space-scanner_mask.jpg']));
    fprintf(fid,'\n\t</div>');
end
fprintf(fid,'\n</div>');
fprintf(fid,'\n</body>');
fprintf(fid,'\n</html>');
fclose(fid);