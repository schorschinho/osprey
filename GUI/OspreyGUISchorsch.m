function OspreyGUI(MRSCont)
%% OspreyGUI(MRSCont)
%   This function creates a one-in-all figure with visualizations of the
%   processed data (spectra in the frequency domain), voxel coregistration
%   and segmentation, and quantification tables.
%
%   The figure contains several tabs, not all of which may be available at
%   all times:
%       - Data
%       - Coregistration and segmentation
%       - Fit
%       - Quantification
%
%   As an example, if coregistration and segmentation have not been
%   performed, the respective tab will be grayed out.
%
%   USAGE:
%       OspreyGUI(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-30)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-30: First version of the code.

% Close any remaining open figures
close all;

% Create the user interface for the application and return a
% structure of handles for global use.
gui = struct();
% Create the overall figure
gui.window = figure(...
    'Name', 'Osprey', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'HandleVisibility', 'off', ...
    'Renderer', 'painters');
% Resize such that width is 1.2941 * height (1.2941 is the ratio
% between width and height of standard US letter size (11x8.5 in).
screenSize      = get(0,'ScreenSize');
canvasSize      = screenSize;
canvasSize(4)   = screenSize(4) * 0.9;
canvasSize(3)   = canvasSize(4) * (11/8.5);
canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
set(gui.window, 'Position', canvasSize);

% Create the main horizontal box division between menu (left) and display
% panel tabs (right)
mainLayout = uix.HBox(...
    'Parent', gui.window);
    % Create the left-side menu
    leftMenu = uix.VBox(...
        'Parent',mainLayout, ...
        'Padding',5);
        % Divide into the upper panel containing the Gannet logo button and the
        % lower panel containing the menu buttons
        p1 = uix.Panel(...
            'Parent', leftMenu, ...
            'Padding', 5);
            % Place Gannet logo on a push button
            b_about = uicontrol(p1,'Style','PushButton');
            set(b_about,'Units','Normalized','Position',[0 0 1 1],'BackgroundColor',[0.93 0.93 0.93]);
            [img, ~, ~] = imread('osprey.png', 'BackgroundColor', [0.93 0.93 0.93]);
            [img2] = imresize(img, 0.15);
            set(b_about,'CData', img2);
            logoFcn = @()imread('osprey.png', 'BackgroundColor', [0.93 0.93 0.93]);
logoBanner = uiw.utility.loadIcon(logoFcn);

d = uiw.dialog.About(...
    'Name', 'Osprey',...
    'Version','0.0.1',...
    'Date', 'October 6, 2019',...
    'Timeout', 3,...
    'CustomText', 'Osprey is provided by Johns Hopkins University.',...
    'ContactInfo', 'gabamrs@gmail.com',...
    'LogoCData', logoBanner);

        p2 = uix.VButtonBox(...
            'Parent', leftMenu, ...
            'Padding', 5, ...
            'Spacing', 10, ...
            'ButtonSize', [180 50]);
        set(leftMenu, 'Heights', [-0.2 -0.8]);
        % Data input button
        b_input = uicontrol('Parent', p2,'Style','PushButton','String','Data input');
        set(b_input,'Units','Normalized','Position',[0.1 0.9 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_input,'Callback',@gui_input);
        % Load button
        b_load = uicontrol('Parent', p2,'Style','PushButton','String','Load data','Enable','off');
        set(b_load,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Fit button
        b_fit = uicontrol('Parent', p2,'Style','PushButton','String','Fit data','Enable','off');
        set(b_fit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Coregister button
        b_coreg = uicontrol('Parent', p2,'Style','PushButton','String','CoRegister','Enable','off');
        set(b_coreg,'Units','Normalized','Position',[0.1 0.59 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Segment button
        b_segm = uicontrol('Parent', p2,'Style','PushButton','String','Segment','Enable','off');
        set(b_segm,'Units','Normalized','Position',[0.1 0.51 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Quantify button
        b_quant = uicontrol('Parent', p2,'Style','PushButton','String','Quantify','Enable','off');
        set(b_quant,'Units','Normalized','Position',[0.1 0.43 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % CoRegStandAlone button
        b_coregSA = uicontrol('Parent', p2,'Style','PushButton','String','CoRegStandAlone','Enable','off');
        set(b_coregSA,'Units','Normalized','Position',[0.1 0.28 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % DeIdentify button
        b_deid = uicontrol('Parent', p2,'Style','PushButton','String','DeIdentify','Enable','off');
        set(b_deid,'Units','Normalized','Position',[0.1 0.2 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_deid,'Callback',@gui_deid);
        % Settings button
        b_settings = uicontrol('Parent', p2,'Style','PushButton','String','Settings');
        set(b_settings,'Units','Normalized','Position',[0.1 0.05 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_settings,'Callback',@gui_settings);
        % Exit button
        b_exit = uicontrol('Parent', p2,'Style','PushButton','String','Exit');
        set(b_exit,'Units','Normalized','Position',[0.1 0 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_exit,'Callback',@onExit);
    
    % Create the display panel tab row
    tabs = uix.TabPanel(...
        'Parent', mainLayout, ...
        'Padding', 5, ...
        'FontName', 'Century Gothic', ...
        'FontSize', 16);
    dataTab     = uix.VBox('Parent', tabs, 'Padding', 5);
    coregTab    = uix.Grid('Parent', tabs, 'Spacing', 5);
    fitTab      = uix.VBox('Parent', tabs, 'Padding', 5);
    quantifyTab = uix.Panel('Parent', tabs, 'Padding', 5);
    tabs.TabTitles  = {'Data', 'Coreg / Segment', 'Fit', 'Quantification'};
    tabs.TabWidth   = 200;
    tabs.Selection  = 1;
    tabs.TabEnables = {'off', 'off', 'off', 'off'};
    
    % Now enable the display tabs depending on which processing steps have
    % been completed:
    if MRSCont.flags.didLoadData % Was data loaded at all that can be looked at?
        tabs.TabEnables{1} = 'on';
        
        %%% DATA DISPLAY PANEL %%%
        % Fill the data display panel
        dataStats = uix.Panel(...
            'Parent', dataTab, ...
            'Padding', 5);
        dataPlot = uix.HBox(...
            'Parent', dataTab, ...
            'Padding', 5);
        set(dataTab, 'Heights', [-0.1 -0.9]);
        
        % Create the axes and control slider
        gui.dataAxes = axes('Parent', dataPlot, 'Linewidth', 1);
        set(gui.dataAxes,'TickDir','out','YTick',[],'YTickLabel',[]);
        gui.dataAxes.XRuler.MinorTick = 'on';
        gui.dataAxes.XGrid = 'on';
        gui.dataAxes.XMinorGrid = 'on';
        gui.dataAxes.GridAlpha = 0.5;
        gui.dataAxes.MinorGridAlpha = 0.5;
        gui.dataAxes.XLabel.String = 'Chemical Shift (ppm)';
        P = getPlotDetails(MRSCont);

        
        

        gui.dataControls = uix.HBox(...
            'Parent', dataPlot, ...
            'Padding', 5);
        gui.dataButtonList = uix.VButtonBox(...
            'Parent', gui.dataControls, ...
            'Padding', 5);

        gui.editbox_plotData_lowPPM     = uicontrol('Parent',gui.dataButtonList, 'String', 'Lower PPM', 'Style','edit',...
            'value', 0, 'FontSize', 12, 'FontName', 'Century Gothic');
        gui.editbox_plotData_highPPM    = uicontrol('Parent',gui.dataButtonList, 'String', 'Upper PPM', 'Style','edit',...
            'value', 5, 'FontSize', 12, 'FontName', 'Century Gothic');
        set(gui.dataButtonList, 'ButtonSize', [100 15], 'Spacing', 10);
        
        gui.sl_plotData = uicontrol('Parent',gui.dataControls,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'Units', 'Normalized', 'Position', [0 0 1 1], ...
            'sliderstep',[1 1],...
            'callback',{@sl_plotData_Call,gui});

        if MRSCont.nDatasets == 1
            gui.sl_plotData.Visible = 'off';
        end
        set(dataPlot, 'Widths', [-1 -0.2]);
        %%% /DATA DISPLAY PANEL %%%
    
    end
    if MRSCont.flags.didCoreg % Have coreg/segment masks been created?
        tabs.TabEnables{2} = 'on';
        
        %%% COREG+SEGMENT DISPLAY PANEL %%%
        % Fill the coreg display panel
        coregPlotTra = uix.Panel(...
            'Parent', coregTab, ...
            'Padding', 5);
        coregPlotCor = uix.Panel(...
            'Parent', coregTab, ...
            'Padding', 5);
        coregPlotSag = uix.Panel(...
            'Parent', coregTab, ...
            'Padding', 5);
        coregStats = uix.Panel(...
            'Parent', coregTab, ...
            'Padding', 5);
        
        % Create the axes and control slider
        [img_tra, img_cor, img_sag] = update_plotCoreg(1);
        S.coregAxesTra = axes('Parent', coregPlotTra, 'ActivePositionProperty', 'Position', 'XTickLabel', [], 'YTickLabel', []);
        coregTraImg = imagesc(img_tra, 'Parent', S.coregAxesTra);
        colormap(S.coregAxesTra, 'gray');
        S.coregAxesCor = axes('Parent', coregPlotCor, 'ActivePositionProperty', 'Position', 'XTickLabel', [], 'YTickLabel', []);
        coregCorImg = imagesc(img_cor, 'Parent', S.coregAxesCor);
        colormap(S.coregAxesCor, 'gray');
        S.coregAxesSag = axes('Parent', coregPlotSag, 'ActivePositionProperty', 'Position', 'XTickLabel', [], 'YTickLabel', []);
        coregSagImg = imagesc(img_sag, 'Parent', S.coregAxesSag);
        colormap(S.coregAxesSag, 'gray');
        S.sl_plotCoreg = uicontrol('Parent',coregTab,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'sliderstep',[1 10],...
            'callback',{@sl_plotCoreg_Call,S});
        if MRSCont.nDatasets == 1
            S.sl_plotCoreg.Visible = 'off';
        end
        set(coregTab, 'Widths',[-1 -1 10], 'Heights', [-1 -1]);
        %%% /COREG+SEGMENT DISPLAY PANEL %%%
    
    end
    if MRSCont.flags.didFit % Has data fitting been run?
        
        %%% FIT DISPLAY PANEL %%%
        tabs.TabEnables{3} = 'on';
        % Fill the fit display panel
        fitStats = uix.Panel(...
            'Parent', fitTab, ...
            'Padding', 5);
        fitPlot = uix.HBox(...
            'Parent', fitTab, ...
            'Padding', 5);
        set(fitTab, 'Heights', [-0.1 -0.9]);
        
        % Create the axes
        S.fitAxes = axes('Parent', fitPlot, 'Linewidth', 1);
        S.plotRangePPM = [0 5];
        line = getLineDetails(1, S);
        S.LN = plot(S.fitAxes,MRSCont.processed.A{1}.ppm,real(MRSCont.processed.A{1}.specs),'k','linewidth',1);
        set(S.fitAxes,'xdir','reverse','xlim',S.plotRangePPM,'ylim',[line.data_min line.data_max]);
        set(S.fitAxes,'TickDir','out','YTick',[],'YTickLabel',[]);
        S.fitAxes.XRuler.MinorTick = 'on';
        S.fitAxes.XRuler.MinorTickValues = S.plotRangePPM(1):0.1:S.plotRangePPM(end);
        S.fitAxes.XRuler.TickValues = (S.plotRangePPM(1):0.2:S.plotRangePPM(end));
        S.fitAxes.XGrid = 'on';
        S.fitAxes.XMinorGrid = 'on';
        S.fitAxes.GridAlpha = 0.5;
        S.fitAxes.MinorGridAlpha = 0.5;
        S.fitAxes.XLabel.String = 'Chemical Shift (ppm)';
        
        % Create the controls
        S.fitControls = uix.HBox(...
            'Parent', fitPlot, ...
            'Padding', 5);
        % Create buttons to switch individual plots on/off
        S.fitButtonList = uix.VButtonBox(...
            'Parent', S.fitControls, ...
            'Padding', 5);
        for kk = 1:length(MRSCont.fit.basisSet.name)
            S.fitMetabButton{kk} = uicontrol( 'Parent', S.fitButtonList, 'String', MRSCont.fit.basisSet.name{kk}, 'Style', 'checkbox', 'FontSize', 12, 'FontName', 'Century Gothic', 'Value', 1);
        end
        set(S.fitButtonList, 'ButtonSize', [100 15], 'Spacing', 10);
        
        % Slider to select datasets on the far right
        S.sl_plotFit = uicontrol('Parent',S.fitControls,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'Units', 'Normalized', 'Position', [0 0 1 1], ...
            'sliderstep',[1 1],...
            'callback',{@sl_plotFit_Call,S});
        set(S.fitControls, 'Widths', [-1 10]);
        
        if MRSCont.nDatasets == 1
            S.sl_plotFit.Visible = 'off';
        end
        
        set(fitPlot, 'Widths', [-1 -0.2]);
        %%% /FIT DISPLAY PANEL %%%
        
        
        tabs.TabEnables{4} = 'on';
        
        
    end
    
set(mainLayout, 'Widths', [-0.2 -0.8]);
    

    
    
%-------------------------------------------------------------------------%
    %------- CALLBACK FUNCTIONS FOR THE MAIN MENU -------% 
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.window );
    end % onExit
    %------- /CALLBACK FUNCTIONS FOR THE MAIN MENU ------% 
    
    
    %------- CALLBACK FUNCTIONS FOR THE DATA DISPLAY MENU -------% 
    function P = getPlotDetails(MRSCont)
        % write this one so that P contains all ppm ranges, all plots, and
        % all details (like SNR) that are needed on the data tab
        % get them FOR ALL THE idx so that it can be easily changed by the
        % callbacks
        for idx_p = 1:MRSCont.nDatasets
            whichFields      = fieldnames(MRSCont.processed);
            for idx_p2 = 1:length(whichFields)
                P.(whichFields{idx_p2}).ppm     = MRSCont.processed.(whichFields{idx_p2}){idx_p}.ppm;
                P.(whichFields{idx_p2}).specs   = MRSCont.processed.(whichFields{idx_p2}){idx_p}.specs;
            end
        end
    end


    function [] = sl_plotData_Call(varargin)
        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        idx = get(h,'value');
        
        data_to_plot         = MRSCont.processed.A{idx}.specs;
        line.data_max        = max(real(data_to_plot));
        line.data_min        = min(real(data_to_plot));
        
        line = getLineDetails(idx,S);
                S.plotRangePPM = [0 5];
        line = getLineDetails(1, S);
        S.LN = plot(gui.dataAxes,MRSCont.raw{1}.ppm,real(MRSCont.raw{1}.specs),'k','linewidth',1);
        set(S.dataAxes,'xdir','reverse','xlim',S.plotRangePPM,'ylim',[line.data_min line.data_max]);
                gui.dataAxes.XRuler.MinorTickValues = S.plotRangePPM(1):0.1:S.plotRangePPM(end);
        gui.dataAxes.XRuler.TickValues = (S.plotRangePPM(1):0.2:S.plotRangePPM(end));
        set(S.LN,'xdata', MRSCont.raw{idx}.ppm,'ydata', real(MRSCont.processed.A{idx}.specs))
        set(S.dataAxes,'ylim', [line.data_min line.data_max])
    end
    %------- /CALLBACK FUNCTIONS FOR THE DATA DISPLAY MENU ------%
    
    
    function [] = sl_plotCoreg_Call(varargin)
        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        idx = get(h,'value');
        [img_tra, img_cor, img_sag] = update_plotCoreg(idx);
        coregTraImg = imagesc(img_tra, 'Parent', S.coregAxesTra);
        coregCorImg = imagesc(img_cor, 'Parent', S.coregAxesCor);
        coregSagImg = imagesc(img_sag, 'Parent', S.coregAxesSag);
    end

    function [img_tra, img_cor, img_sag] = update_plotCoreg(idx)
        vol_image   = MRSCont.coreg.vol_image{idx};
        vol_mask    = MRSCont.coreg.vol_mask{idx};
        vox_offs    = [-MRSCont.processed.A{idx}.geometry.pos.PosSag -MRSCont.processed.A{idx}.geometry.pos.PosCor MRSCont.processed.A{idx}.geometry.pos.PosTra];
        T1_max      = MRSCont.coreg.T1_max{idx};
        [img_tra, img_cor, img_sag]       = voxel2world_space(vol_image, vox_offs);
        [mask_tra, mask_cor, mask_sag]    = voxel2world_space(vol_mask, vox_offs);
        
        img_tra = flipud(img_tra/T1_max);
        img_cor = flipud(img_cor/T1_max);
        img_sag = flipud(img_sag/T1_max);
        
        img_tra = img_tra + 0.25*flipud(mask_tra);
        img_cor = img_cor + 0.25*flipud(mask_cor);
        img_sag = img_sag + 0.25*flipud(mask_sag);

    end

    function [] = sl_plotFit_Call(varargin)
        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        idx = get(h,'value');
        line = getLineDetails(idx,S);
        set(S.LN,'xdata', MRSCont.raw{idx}.ppm,'ydata', real(MRSCont.processed.A{idx}.specs))
        set(S.fitAxes,'ylim', [line.data_min line.data_max])
    end

    function [] = extractFit(MRSCont, kk, ll)
        % Initialise variables
        basisSet            = MRSCont.fit.basisSet;
        nMets               = basisSet.nMets;
        nMM                 = basisSet.nMM;
        nBasisFcts          = nMets + nMM;
        N_b                 = length(MRSCont.fit.output.knotLocations);
        fitRangePPM         = MRSCont.fit.settings.fitRange;
        
        % Prepare data, basis sets, result field names, and spline knot spacing to
        % be extracted
        if MRSCont.flags.isUnEdited
            dataFields      = {'A'};
            data = op_ampScale(MRSCont.processed.A{kk}, 1/MRSCont.fit.scale{kk});
            resultFields    = {'resultsUnEdited'};
        elseif MRSCont.flags.isMEGA
            dataFields      = {'diff1','A'};
            data = [real(MRSCont.fit.diff1{kk}.specs); real(MRSCont.fit.A{kk}.specs); imag(MRSCont.fit.diff1{kk}.specs); imag(MRSCont.fit.A{kk}.specs)];
            switch MRSCont.fit.settings.fitStyle
                case 'Concatenated'
                    resultFields = {'resultsMEGA'};
                case 'Separate'
                    resultFields = {'resultsDiff', 'resultsOff'};
            end
        elseif MRSCont.flags.isHERMES
            dataFields      = {'diff1', 'diff2', 'A'};
            data = [real(MRSCont.fit.diff1{kk}.specs); real(MRSCont.fit.diff2{kk}.specs); real(MRSCont.fit.A{kk}.specs); imag(MRSCont.fit.diff1{kk}.specs); imag(MRSCont.fit.diff2{kk}.specs); imag(MRSCont.fit.A{kk}.specs)];
            switch MRSCont.fit.settings.fitStyle
                case 'Concatenated'
                    resultFields = {'resultsHERMES'};
                case 'Separate'
                    resultFields = {'resultsDiff1', 'resultsDiff2', 'resultsOff'};
            end
        elseif MRSCont.flags.isHERCULES
            dataFields      = {'diff1', 'diff2', 'sum'};
            data = [real(MRSCont.fit.diff1{kk}.specs); real(MRSCont.fit.diff2{kk}.specs); real(MRSCont.fit.sum{kk}.specs); imag(MRSCont.fit.diff1{kk}.specs); imag(MRSCont.fit.diff2{kk}.specs); imag(MRSCont.fit.sum{kk}.specs)];
            switch MRSCont.fit.settings.fitStyle
                case 'Concatenated'
                    resultFields = {'resultsHERC'};
                case 'Separate'
                    resultFields = {'resultsDiff1', 'resultsDiff2', 'resultsSum'};
            end
        end
        
        % Select the correct field to extract the fit parameters from,
        % depending on whether the fit was performed on all sub-spectra
        % simultaneously or separately.
        switch MRSCont.fit.settings.fitStyle
            case 'Concatenated'
                ampl        = MRSCont.fit.(resultFields{1}).ampl{kk}; % amplitudes
                zeroPhase   = MRSCont.fit.(resultFields{1}).ph0{kk}; % zero-order phase correction [rad]
                firstPhase  = MRSCont.fit.(resultFields{1}).ph1{kk}; % first-order phase correction [rad]
                gaussLB     = MRSCont.fit.(resultFields{1}).gaussLB{kk}; % Gaussian damping [Hz^2]
                lorentzLB   = MRSCont.fit.(resultFields{1}).lorentzLB{kk}; % Lorentzian damping [Hz] for each basis function
                freqShift   = MRSCont.fit.(resultFields{1}).freqShift{kk}; % Frequency shift [Hz] for each basis function
                spl_pos     = MRSCont.fit.(resultFields{1}).spl_pos{kk}; % spline positions (pts)
                beta_j      = MRSCont.fit.(resultFields{1}).beta_j{kk}; % spline coefficients
                tCrRatio    = MRSCont.fit.(resultFields{1}).tCrRatio{kk}; % tCr ratios
            case 'Separate'
                ampl        = MRSCont.fit.(resultFields{ll}).ampl{kk}; % amplitudes
                zeroPhase   = MRSCont.fit.(resultFields{ll}).ph0{kk}; % zero-order phase correction [rad]
                firstPhase  = MRSCont.fit.(resultFields{ll}).ph1{kk}; % first-order phase correction [rad]
                gaussLB     = MRSCont.fit.(resultFields{ll}).gaussLB{kk}; % Gaussian damping [Hz^2]
                lorentzLB   = MRSCont.fit.(resultFields{ll}).lorentzLB{kk}; % Lorentzian damping [Hz] for each basis function
                freqShift   = MRSCont.fit.(resultFields{ll}).freqShift{kk}; % Frequency shift [Hz] for each basis function
                spl_pos     = MRSCont.fit.(resultFields{ll}).spl_pos{kk}; % spline positions (pts)
                beta_j      = MRSCont.fit.(resultFields{ll}).beta_j{kk}; % spline coefficients
                tCrRatio    = MRSCont.fit.(resultFields{ll}).tCrRatio{kk}; % tCr ratios
        end
        
        % 3. Run the time-domain operations on the basis functions
        % (frequency shift, Lorentzian damping, Gaussian damping)
        t = basisSet.t;
        for ii=1:nBasisFcts
            basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp([1i*freqShift(ii) - lorentzLB(ii) - gaussLB.*t].*t)';
        end
        basisSet.specs = fftshift(ifft(basisSet.fids,[],1),1);
        
        % 4. Run the frequency-domain operations on the basis functions
        % (zero and first order phase correction)
        f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
        for ii=1:nBasisFcts
            basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp([1i*zeroPhase + 1i*firstPhase*2*pi*f])';
        end
        basisSet.fids = fft(fftshift(basisSet.specs,1),[],1);
        % Cut out the frequency range of the basis set
        basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end));
        
        % 5. Set up baseline spline
        test = spmak(spl_pos, beta_j');
        B = fnval(test,1:1:basisSet.sz(1));
        % % Add baseline only to the real part of the spectrum
        % B = [B zeros(1,length(B))]';
        % Get imaginary part through Hilbert transform
        B_Hilb = hilbert(B);
        B = [real(B_Hilb)+1i*imag(B_Hilb)]';
        
        % 6. Cut out the frequency range of the spectrum to be fit
        dataToFit   = op_freqrange(data,fitRangePPM(1),fitRangePPM(end));
        A           = basisSet.specs;
        fitted      = (A*ampl + B);
        resid       = dataToFit.specs - fitted;
        
        % Select the data range to plot
        data_to_plot    = dataToFit.specs;
        fitted_to_plot  = fitted;
        resid_to_plot   = resid;
        B_to_plot       = B;

    end
end

