classdef OspreyGUIapp < handle
%% OspreyGUIapp
%   This class creates a one-in-all figure with visualizations of the
%   processed data (spectra in the frequency domain), voxel coregistration
%   and segmentation, quantification tables, and results.
%
%   The figure contains several tabs, not all of which may be available at
%   all times:
%       - Data
%       - Processed
%       - Coregistration and segmentation
%       - Fit
%       - Quantification
%       - Overview
%
%   As an example, if coregistration and segmentation have not been
%   performed, the respective tab will be grayed out.
%
%   USAGE:
%       OspreyGUIapp;
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2019-11-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-13: First version of the code.

 %% Properties
 % Here the properties for the gui class are defined
    properties
        MRSCont % MRS container with data
        folder % all folders
        colormap % colormaps for gui and 
        controls % gui control handles
        load % OspreyLoad infos
        process %  OpsreyProcess infos
        fit % OpsreyFit infos
        quant % OspreyQuant infos
        overview % OspreyOverview infos
        figure % figure handle
        layout % Layout infos
        upperBox % contains upper part of gui
        Plot % plot struct
        InfoText % info text struct (left side of fit plot)
        Results % result struct (fit results)
    end
    
    methods
        
        function gui = OspreyGUIapp(MRSCont)
        % This is the "constructor" for the class
        % It runs when an object of this class is created
        %% Initialize variables
        % Close any remaining open figures & add folders
            close all;
            [settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
            gui.folder.allFolders      = strsplit(settingsFolder, filesep);
            gui.folder.ospFolder       = strjoin(gui.folder.allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
        % SPM
            if isfile(fullfile(gui.folder.ospFolder,'GUI','SPMpath.mat'))
                load(fullfile(gui.folder.ospFolder,'GUI','SPMpath.mat'),'SPMpath')
                gui.folder.spmversion = SPMpath;
            else
                gui.folder.spmversion = uipickfiles('FilterSpec',[gui.folder.ospFolder],'REFilter', '\$','NumFiles',1,'Prompt','Select your SPM-folder (Will be saved in SPMpath.mat file in the GUI folder)');
                gui.folder.spmversion = gui.folder.spmversion{1};
                SPMpath = gui.folder.spmversion;
                save(fullfile(gui.folder.ospFolder,'GUI','SPMpath.mat'),'SPMpath');
            end
        % Check if SPM12 is installed
            if isempty(gui.folder.spmversion)
                error('SPM not found! Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.');
            elseif strcmpi(gui.folder.spmversion(end-3:end),'spm8')
                error(['SPM8 detected, but only SPM12 is supported. ' ...
                       'Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.']);
            end
       % Load selected colormap     
            gui.colormap = MRSCont.colormap;
            
        %Setting up inital values for the gui control variables
        %Global controls
            gui.controls.Selected = 1;
            gui.controls.Number = 1;
            gui.controls.KeyPress = 0;
        %File selections for each sub function        
            gui.load.Selected = 1;
            gui.process.Selected = 1;
            gui.fit.Selected = 1;
            gui.quant.Selected.Model = 1;
            gui.quant.Selected.Quant = 1;
            gui.overview.Selected.Metab = 1;
        %Names for each selection    
            gui.load.Names.Spec = {'metabolites'};
            
        %Setting up remaining values in dependence of the conducted processing steps
            if MRSCont.flags.didLoadData %Get variables regarding the loading
                if ~isempty(MRSCont.raw{1,gui.controls.Selected}.seq)
                    if strcmp(sprintf('\n'),MRSCont.raw{1,gui.controls.Selected}.seq(end)) %Clean up Sequence Name if needed
                        gui.load.Names.Seq = MRSCont.raw{1,gui.controls.Selected}.seq(1:end-1);
                    else
                        gui.load.Names.Seq = MRSCont.raw{1,gui.controls.Selected}.seq;
                    end
                else
                    if MRSCont.flags.isUnEdited
                        gui.load.Names.Seq =['Unedited ' MRSCont.vendor];
                    end
                    if MRSCont.flags.isMEGA
                        gui.load.Names.Seq =['MEGA ' MRSCont.vendor];
                    end
                    if MRSCont.flags.isHERMES
                        gui.load.Names.Seq =['HERMES ' MRSCont.vendor];
                    end
                    if MRSCont.flags.isHERCULES
                        gui.load.Names.Seq =['HERCULES ' MRSCont.vendor];
                    end
                    if MRSCont.flags.isPRIAM
                        gui.load.Names.Seq =['PRIAM ' MRSCont.vendor];
                    end
                end
                gui.load.Names.Geom = fieldnames(MRSCont.raw{1,1}.geometry.size); %Get variables regarding voxel geometry
            end
            
            if MRSCont.flags.didProcess %Get variables regarding the processing
                gui.process.Number = length(fieldnames(MRSCont.processed));
                gui.process.Names = fieldnames(MRSCont.processed);
                gui.process.Names = sort(gui.process.Names);
            end         

            if MRSCont.flags.didFit %Get variables regarding the fitting
                if strcmp(MRSCont.opts.fit.style, 'Concatenated')
                    temp = fieldnames(MRSCont.fit.results);
                    if MRSCont.flags.isUnEdited
                        gui.fit.Names = fieldnames(MRSCont.fit.results);
                    end
                    if MRSCont.flags.isMEGA
                        gui.fit.Names = {'diff1','sum'};
                        if length(temp) == 2
                            gui.fit.Names{3} = temp{2};
                        else if length(temp) == 3
                            gui.fit.Names{3} = temp{2};
                            gui.fit.Names{4} = temp{3};
                            end
                        end
                    end
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                    gui.fit.Names = {'diff1','diff2','sum'};
                    if length(temp) == 2
                        gui.fit.Names{4} = temp{2};
                    else if length(temp) == 3
                        gui.fit.Names{4} = temp{2};
                        gui.fit.Names{5} = temp{3};
                        end
                    end
                end
                    gui.fit.Number = length(gui.fit.Names); 
                else
                    gui.fit.Names = fieldnames(MRSCont.fit.results);
                    gui.fit.Number = length(fieldnames(MRSCont.fit.results));   
                end
            end

            if MRSCont.flags.didQuantify %Get variables regarding the quantification
                gui.quant.Number.Model = length(fieldnames(MRSCont.quantify.tables));
                gui.quant.Names.Model = fieldnames(MRSCont.quantify.tables);
                gui.quant.Number.Quants = length(fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1})));
                gui.quant.Names.Quants = fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1}));
                gui.quant.Number.Metabs = length(MRSCont.quantify.metabs);
                gui.quant.Selected.Metab = find(strcmp(MRSCont.quantify.metabs, 'tNAA'));
                gui.quant.Selected.Model = 1;
                gui.quant.idx.GABA = find(strcmp(MRSCont.quantify.metabs, 'GABA'));
            end
            if MRSCont.flags.didOverview %Get variables for the overview tab
                gui.overview.NAAnormed = 1;
                gui.overview.Number.Groups = MRSCont.overview.NoGroups;
                [gui.colormap.cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
                temp = gui.colormap.cb(3,:);
                gui.colormap.cb(3,:) = gui.colormap.cb(4,:);
                gui.colormap.cb(4,:) = temp;
                if isfield(MRSCont.overview, 'corr')
                    gui.overview.Names.Corr = MRSCont.overview.corr.Names{1};
                    gui.overview.CorrMeas = MRSCont.overview.corr.Meas;
                end
                gui.overview.Selected.Corr = 1;
                gui.overview.Selected.CorrChoice = 1;
                gui.overview.Names.QM = {'SNR','FWHM (ppm)'};
            end
            
        %% Create the overall figure
            gui.figure = figure('Name', 'Osprey', 'NumberTitle', 'off', 'Visible', 'on',...
                                'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', gui.colormap.Background);
            setappdata(gui.figure,'MRSCont',MRSCont);
        % Resize such that width is 1.2941 * height (1.2941 is the ratio
        % between width and height of standard US letter size (11x8.5 in).
            screenSize      = get(0,'ScreenSize');
            canvasSize      = screenSize;
            canvasSize(4)   = screenSize(4) * 0.9;
            canvasSize(3)   = canvasSize(4) * (11/8.5);
            canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
            canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
            set(gui.figure, 'Position', canvasSize);

        % Create the main horizontal box division between menu (left) and display
        % panel tabs (right)
            gui.layout.mainLayout = uix.HBox('Parent', gui.figure,'BackgroundColor',gui.colormap.Background);
        % Create the left-side menu
            gui.layout.leftMenu = uix.VBox('Parent',gui.layout.mainLayout, 'Padding',4,'Spacing', 2,'BackgroundColor',gui.colormap.Background);
            gui.layout.Buttonbox = uix.HBox('Parent',gui.layout.leftMenu, 'BackgroundColor',gui.colormap.Background);
            gui.layout.b_about = uicontrol(gui.layout.Buttonbox,'Style','PushButton');
            [img, ~, ~] = imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
            [img2] = imresize(img, 0.08);
            set(gui.layout.b_about,'CData', img2, 'TooltipString', 'Contact developers via mail');
            set(gui.layout.b_about,'Callback',{@osp_onOsp});
           
           
           gui.layout.Infobuttons = uix.VButtonBox('Parent',gui.layout.Buttonbox,'BackgroundColor',gui.colormap.Background);
        % Divide into the upper panel containing the Osprey logo button

           
           gui.layout.b_pub = uicontrol(gui.layout.Infobuttons,'Style','PushButton');
           [img, ~, ~] = imread('PubMed.png', 'BackgroundColor', gui.colormap.Background);
           [img2] = imresize(img, 0.06);
           set(gui.layout.b_pub,'CData', img2, 'TooltipString', 'Cite these papers when using Osprey');
           set(gui.layout.b_pub,'Callback',{@osp_onPub});
           
            gui.layout.b_git = uicontrol(gui.layout.Infobuttons,'Style','PushButton');
           [img, ~, ~] = imread('GitHubB.png', 'BackgroundColor', gui.colormap.Background);
           [img2] = imresize(img, 0.07);
           set(gui.layout.b_git,'CData', img2, 'TooltipString', 'Keep yourself updated and request/develop new features on GitHub');
           set(gui.layout.b_git,'Callback',{@osp_onGit});

           set(gui.layout.Infobuttons, 'ButtonSize', [150 100]);
           set(gui.layout.Buttonbox, 'Width', [-0.5 -0.5]);     
        %% Create left menu
            gui.layout.p2 = uix.VButtonBox('Parent', gui.layout.leftMenu, 'Spacing', 3, ...
                            'BackgroundColor',gui.colormap.Background);
            set(gui.layout.leftMenu, 'Heights', [-0.2 -0.8]);
            set(gui.layout.p2, 'ButtonSize', [300 60]);
        % Load button
            gui.layout.b_load = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Load data','Enable','on','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_load,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if MRSCont.flags.didLoadData
                gui.layout.b_load.Enable = 'off';
            end
            set(gui.layout.b_load,'Callback',{@osp_onLoad,gui}, 'TooltipString', 'Call OspreyLoad');
        % Process button
            gui.layout.b_proc = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Process data','Enable','on','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_proc,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if MRSCont.flags.didProcess
                gui.layout.b_proc.Enable = 'off';
            else if ~MRSCont.flags.didLoadData
                    gui.layout.b_proc.Enable = 'off';
                end 
            end
            set(gui.layout.b_proc,'Callback',{@osp_onProc,gui}, 'TooltipString', 'Call OspreyProcess');
        % Fit button
            gui.layout.b_fit = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Fit data','Enable','on','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_fit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if MRSCont.flags.didFit
                gui.layout.b_fit.Enable = 'off';
            else if ~MRSCont.flags.didProcess
                    gui.layout.b_fit.Enable = 'off';
                end
            end
            set(gui.layout.b_fit,'Callback',{@osp_onFit,gui}, 'TooltipString', 'Call OspreyFit');
        % Coregister button
            gui.layout.b_coreg = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','CoRegister','Enable','off','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_coreg,'Units','Normalized','Position',[0.1 0.59 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if ~isempty(MRSCont.files_nii) && ~MRSCont.flags.didCoreg && MRSCont.flags.didLoadData
                gui.layout.b_coreg.Enable = 'on';
            end
            set(gui.layout.b_coreg,'Callback',{@osp_onCoreg,gui}, 'TooltipString', 'Call OspreyCoreg');
        % Segment button
            gui.layout.b_segm = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Segment','Enable','off','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_segm,'Units','Normalized','Position',[0.1 0.51 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if ~isempty(MRSCont.files_nii) && ~MRSCont.flags.didSeg && MRSCont.flags.didCoreg
                gui.layout.b_segm.Enable = 'on';
            end
            set(gui.layout.b_segm,'Callback',{@osp_onSeg,gui}, 'TooltipString', 'Call OspreySeg');
        % Quantify button
            gui.layout.b_quant = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Quantify','Enable','on','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_quant,'Units','Normalized','Position',[0.1 0.43 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            if MRSCont.flags.didQuantify
                gui.layout.b_quant.Enable = 'off';
            else if ~MRSCont.flags.didFit
                    gui.layout.b_quant.Enable = 'off';
                end
            end
            set(gui.layout.b_quant,'Callback',{@osp_onQuant,gui}, 'TooltipString', 'Call OspreyQuantify');
        % DeIdentify button
            gui.layout.b_deid = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','DeIdentify','Enable','off','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_deid,'Units','Normalized','Position',[0.1 0.2 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            set(gui.layout.b_deid,'Callback',{@gui_deid,gui,MRSCont}, 'TooltipString', 'DeIndentify');
        % Save button
            gui.layout.b_save = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Save MRSCont','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_save,'Units','Normalized','Position',[0.1 0 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            set(gui.layout.b_save,'Callback',{@osp_onSave,gui}, 'TooltipString', 'Save MRSCont as .mat-file');
        % Exit button
            gui.layout.b_exit = uicontrol('Parent', gui.layout.p2,'Style','PushButton','String','Exit','ForegroundColor', gui.colormap.Foreground);
            set(gui.layout.b_exit,'Units','Normalized','Position',[0.1 0 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
            set(gui.layout.b_exit,'Callback',{@osp_onExit,gui}, 'TooltipString', 'See you next time!');
        % Create list of files for the Listbox
            gui.layout.controlPanel = uix.Panel('Parent', gui.layout.leftMenu, 'Title', 'MRS Container','BackgroundColor',gui.colormap.Background);
            set(gui.layout.controlPanel,'Units','Normalized','Position',[0.5 0 0.66 0.1], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold', 'ForegroundColor',gui.colormap.Foreground, 'HighlightColor',gui.colormap.Foreground, 'ShadowColor',gui.colormap.Foreground);
            gui.layout.fileList = MRSCont.files;
            SepFileList = cell(1,length(MRSCont.files));
            gui.layout.RedFileList = cell(1,length(MRSCont.files));
            gui.layout.OnlyFileList = cell(1,length(MRSCont.files));
            for i = 1 : length(MRSCont.files) %find last two subfolders and file names
                SepFileList{i} =  split(gui.layout.fileList(i), filesep);
                if length(SepFileList{i}) == 1
                    SepFileList{i} =  split(gui.layout.fileList(i), '\');
                end
                gui.layout.RedFileList{i} = [filesep SepFileList{i}{end-2} filesep SepFileList{i}{end-1} filesep SepFileList{i}{end}];
                gui.layout.OnlyFileList{i} = [SepFileList{i}{end}];
            end
            clear SepFileList
            gui.layout.ListBox = uicontrol('Style', 'list','BackgroundColor', 'w','FontName', 'Arial','BackgroundColor',gui.colormap.Background, ...
                                    'Parent', gui.layout.controlPanel, 'String',gui.layout.RedFileList(:) , ...
                                    'Value', gui.controls.Selected, 'Interruptible', 'on', 'BusyAction', 'cancel', ...
                                    'ForegroundColor',gui.colormap.Foreground, 'TooltipString', 'Select a file you want to inspect.');
        %% Create the display panel tab row

            gui.layout.tabs = uix.TabPanel('Parent', gui.layout.mainLayout, 'Padding', 3, 'FontName', 'Arial',...
                            'FontSize', 16,'BackgroundColor',gui.colormap.Background,...
                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.layout.rawTab      = uix.TabPanel('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background,...
                                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                            'FontName', 'Arial', 'TabLocation','bottom','FontSize', 10);
            gui.layout.proTab      = uix.TabPanel('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background,...
                                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                        'FontName', 'Arial', 'TabLocation','bottom','FontSize', 10);
            gui.layout.fitTab      = uix.TabPanel('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background,...
                                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                        'FontName', 'Arial', 'TabLocation','bottom','FontSize', 10);                        
            gui.layout.coregTab    = uix.VBox('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background);
            gui.layout.quantifyTab = uix.TabPanel('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background,...
                                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                        'FontName', 'Arial', 'TabLocation','bottom','FontSize', 10);
            gui.layout.overviewTab = uix.TabPanel('Parent', gui.layout.tabs, 'BackgroundColor',gui.colormap.Background,...
                                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                        'FontName', 'Arial', 'TabLocation','bottom','FontSize', 10);
            gui.layout.tabs.TabTitles  = {'Raw', 'Processed', 'Fit', 'Cor/Seg', 'Quantified','Overview'};
            gui.layout.tabs.TabWidth   = 115;
            gui.layout.tabs.Selection  = 1;
            gui.layout.tabs.TabEnables = {'off', 'off', 'off', 'off', 'off', 'off'};
            set( gui.layout.mainLayout, 'Widths', [-0.2  -0.8] );
            gui.layout.EmptydataPlot = 0;
        %% Here we create the inital setup of the tabs
        % Now enable the display tabs depending on which processing steps have
        % been completed:
            if MRSCont.flags.didLoadData % Was data loaded at all that can be looked at?
                osp_iniLoadWindow(gui);
                set(gui.controls.b_save_RawTab,'Callback',{@osp_onPrint,gui});
            end
            if MRSCont.flags.didProcess % Has data fitting been run?
                osp_iniProcessWindow(gui);
                set(gui.controls.b_save_proTab,'Callback',{@osp_onPrint,gui});
            end
            if MRSCont.flags.didFit % Has data fitting been run?
                osp_iniFitWindow(gui);
                set(gui.controls.b_save_fitTab,'Callback',{@osp_onPrint,gui});
            end            
            if MRSCont.flags.didCoreg % Have coreg/segment masks been created?
                osp_iniCoregWindow(gui);
                set(gui.controls.b_save_coregTab,'Callback',{@osp_onPrint,gui});
            end
            if MRSCont.flags.didQuantify % Has data fitting been run?
                osp_iniQuantifyWindow(gui);
            end
            if MRSCont.flags.didOverview % Has data fitting been run?
                osp_iniOverviewWindow(gui);
                set(gui.layout.overviewTab, 'SelectionChangedFcn',{@osp_OverviewTabChangedFcn,gui});
                set(gui.controls.pop_specsOvPlot,'callback',{@osp_pop_specsOvPlot_Call,gui});
                set(gui.controls.pop_meanOvPlot,'callback',{@osp_pop_meanOvPlot_Call,gui});
                set(gui.controls.pop_quantOvPlot,'callback',{@osp_pop_quantOvPlot_Call,gui});
                set(gui.controls.pop_distrOvQuant,'callback',{@osp_pop_distrOvQuant_Call,gui});
                set(gui.controls.pop_distrOvMetab,'callback',{@osp_pop_distrOvMetab_Call,gui});
                set(gui.controls.pop_corrOvQuant,'callback',{@osp_pop_corrOvQuant_Call,gui});
                set(gui.controls.pop_corrOvMetab,'callback',{@osp_pop_corrOvMetab_Call,gui});
                set(gui.controls.pop_corrOvCorr,'callback',{@osp_pop_corrOvCorr_Call,gui});
                set(gui.controls.pop_whichcorrOvCorr,'callback',{@osp_pop_whichcorrOvCorr_Call,gui});
                set(gui.controls.b_save_specOvTab,'Callback',{@osp_onPrint,gui});
                set(gui.controls.b_save_meanOvTab,'Callback',{@osp_onPrint,gui});
                set(gui.controls.b_save_distrOvTab,'Callback',{@osp_onPrint,gui});
                set(gui.controls.b_save_corrOvTab,'Callback',{@osp_onPrint,gui});
            end
            gui.layout.tabs.Selection  = 1;
            if ~MRSCont.flags.didLoadData %Turn of Listbox if data has not been loaded    
                gui.layout.ListBox.Enable = 'off';
            end
            set(gui.figure, 'Visible','on');

        %% Here we add callback listeners triggered on selection changes
            set(gui.layout.tabs,'SelectionChangedFcn',{@osp_SelectionChangedFcn,gui});
            set(gui.layout.rawTab, 'SelectionChangedFcn',{@osp_RawTabChangeFcn,gui});
            set(gui.layout.proTab,'SelectionChangedFcn',{@osp_ProTabChangeFcn,gui});
            set(gui.layout.fitTab, 'SelectionChangedFcn',{@osp_FitTabChangeFcn,gui});
            set(gui.layout.quantifyTab, 'SelectionChangedFcn',{@osp_QuantTabChangeFcn,gui});            
            set(gui.layout.ListBox,'Callback', {@osp_onListSelection,gui},'KeyPressFcn',{@osp_WindowKeyDown,gui}, 'KeyReleaseFcn', {@osp_WindowKeyUp,gui});     
        end                
    end
end                                                      % End of class definition