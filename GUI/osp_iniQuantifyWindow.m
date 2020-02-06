function osp_iniQuantifyWindow(gui)
%% osp_iniQuantifyWindow
%   This function creates the inital quantify window in the gui.
%
%
%   USAGE:
%       osp_iniQuantifyWindow(gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-16: First version of the code.
%%% 1. GET HANDLES %%%

% This function creates the initial quantify window
        MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
        gui.layout.tabs.TabEnables{5} = 'on';
        gui.layout.tabs.Selection  = 5;
        gui.layout.EmptyQuantPlot = 0;
        
        %%% 2. CREATING SUB TABS FOR THIS TAB %%%%
% In this case one tab fo each fit (off,sum,diff1,diff2,ref,water)
         gui.layout.quantifyTab.TabWidth   = 115;
         for t = 1 : gui.quant.Number.Model %Create tabs depending on the number of fits
                gui.layout.(['quantTab' gui.quant.Names.Model{t}]) = uix.VBox('Parent', gui.layout.quantifyTab,...
                                                                               'BackgroundColor',gui.colormap.Background);
                gui.layout.quantifyTabhandles{t} = ['quantTab' gui.quant.Names.Model{t}];
         end
        gui.layout.quantifyTab.TabTitles  = gui.quant.Names.Model;
        
%%% 3. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        for t = 1 : gui.quant.Number.Model %Loop over fits
            gui.upperBox.quant.Info = uix.Panel('Parent', gui.layout.(gui.layout.quantifyTabhandles{t}), 'Padding', 5, ...
                                      'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                                      'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,...
                                      'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            % Creates layout for plotting and data control
            gui.Plot.quant = uix.HBox('Parent', gui.layout.(gui.layout.quantifyTabhandles{t}),'BackgroundColor',gui.colormap.Background);
            set(gui.layout.(gui.layout.quantifyTabhandles{t}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
           StatText = ['Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style ...
                         '\nSelected subspecs: ' gui.quant.Names.Model{gui.quant.Selected.Model} ];
                     gui.InfoText.quant  = uicontrol('Parent',gui.upperBox.quant.Info,'style','text',...
                'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
%%% 4. VISUALIZATION PART OF THIS TAB %%%%
% In this case a table is created based on a uicontol slider
            gui.quant.Number.Quants = length(fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{t})));
            gui.quant.Names.Quants = fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{t}));
            QuantText = cell(length(MRSCont.quantify.metabs)+1,gui.quant.Number.Quants);
            QuantText{1,1} = 'Metabolite';
            QuantText(2:end,1) = MRSCont.quantify.metabs';
                for q = 1 : gui.quant.Number.Quants % Collect all results
                    QuantText(1,q+1) = gui.quant.Names.Quants(q);
                    if (strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaled') || strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaledGroupNormed')) && isfield(MRSCont.quantify.tables.(gui.quant.Names.Model{t}),'AlphaCorrWaterScaled')
                        idx_GABA  = find(strcmp(MRSCont.quantify.metabs,'GABA'));
                        tempQuantText = cell(length(MRSCont.quantify.metabs),1);
                        tempQuantText(idx_GABA) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{t}).(gui.quant.Names.Quants{q})(gui.controls.Selected,:))';
                        QuantText(2:end,q+1) = tempQuantText;
                    else
                        QuantText(2:end,q+1) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{t}).(gui.quant.Names.Quants{q})(gui.controls.Selected,:))';
                    end
                end
            temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantText,...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            set ( temp, 'BackgroundColor',gui.colormap.Background);
            set( temp, 'Parent', gui.Plot.quant);
        end
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end