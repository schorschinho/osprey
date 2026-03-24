function osp_updateQuantifyWindow(gui)
%% osp_updateQuantifyWindow
%   This function updates the quantify tab.
%
%
%   USAGE:
%       osp_updateQuantifyWindow(gui);
%
%   INPUT:  
%           gui      = gui class containing all handles and the MRSCont             
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
%%% 1. INITIALIZE %%%
        MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class
        gui.layout.EmptyQuantPlot = 0;
        if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            gui.upperBox.quant.Info = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2);
            gui.Plot.quant = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(1);
            gui.InfoText.quant = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children;
            StatText = ['Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style ...
                         '\nSelected subspecs: ' gui.quant.Names.Model{gui.quant.Selected.Model} ];
        elseif  (MRSCont.flags.isPRIAM && isfield(MRSCont.flags,'isPRIAM')) 
            gui.upperBox.quant.Info = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(1);          
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(6),'String',gui.controls.act_x)
            gui.Plot.quant = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(1);
            gui.InfoText.quant = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(1).Children;
            StatText = ['Voxel ' num2str(gui.controls.act_x) ': Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style ...
                         '\nSelected subspecs: ' gui.quant.Names.Model{gui.quant.Selected.Model} ];
        else
            gui.upperBox.quant.Info = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(1);          
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(2).Children.Children.Children(6),'String',gui.controls.act_x)
            gui.Plot.quantMainBox = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(1);            
            gui.Plot.quantHBox = gui.Plot.quantMainBox.Children(1);
            gui.Plot.quantMRSImap = gui.Plot.quantHBox.Children(1).Children(2);
            gui.InfoText.quant = gui.layout.(gui.layout.quantifyTabhandles{gui.quant.Selected.Model}).Children(2).Children(1).Children;
            StatText = ['Metabolite maps and fits of ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style ...
                         '\nSelected subspecs: ' gui.quant.Names.Model{gui.quant.Selected.Model} 'Selected voxel for fit display ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y)];
        end
     
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here           
        set(gui.InfoText.quant, 'String',sprintf(StatText))
%%% 3. VISUALIZATION PART OF THIS TAB %%%
    if ~(isfield(MRSCont.flags,'isMRSI')&& MRSCont.flags.isMRSI)
        gui.quant.Number.Quants = length(fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model})));
        gui.quant.Names.Quants = fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}));
        QuantText = cell(length(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model}))+1,gui.quant.Number.Quants);
        QuantText{1,1} = 'Metabolite';
        QuantText(2:end,1) = MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model})';
    end
       if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
           for q = 1 : gui.quant.Number.Quants %Collect all results
                QuantText(1,q+1) = gui.quant.Names.Quants(q);
                if strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaled') || strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaledGroupNormed')
                    idx_GABA  = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model}),'GABA'));
                    if strcmp(MRSCont.opts.fit.coMM3, 'none')                            
                         tempQuantText = cell(length(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model})),1);
                         tempQuantText(idx_GABA) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).Voxel_1(gui.controls.Selected,:))';
                    else                              
                         tempQuantText = cell(length(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model})),1);
                         tempQuants = MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).Voxel_1(gui.controls.Selected,:);
                         tempQuantText(idx_GABA) = table2cell(tempQuants(1,1));
                         idx_GABAp  = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model}),'GABAplus'));
                         tempQuantText(idx_GABAp) = table2cell(tempQuants(1,2));
                    end  
                    QuantText(2:end,q+1) = tempQuantText;
                else
                    QuantText(2:end,q+1) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).Voxel_1(gui.controls.Selected,:))';
                end
            end
            temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantText,...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
             set( temp, 'BackgroundColor',gui.colormap.Background);
            delete(gui.Plot.quant.Children)
            set( temp, 'Parent', gui.Plot.quant ); %Update table
            set(gui.upperBox.quant.Info,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}]); %Update info Title
        elseif  (MRSCont.flags.isPRIAM && isfield(MRSCont.flags,'isPRIAM')) 
              for q = 1 : gui.quant.Number.Quants %Collect all results
                QuantText(1,q+1) = gui.quant.Names.Quants(q);
                if strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaled') || strcmp(gui.quant.Names.Quants(q),'AlphaCorrWaterScaledGroupNormed')
                    idx_GABA  = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model}),'GABA'));
                    if strcmp(MRSCont.opts.fit.coMM3, 'none')                            
                         tempQuantText = cell(length(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model})),1);
                         tempQuantText(idx_GABA) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).(['Voxel_' num2str(gui.controls.act_x)])(gui.controls.Selected,:))';
                    else                              
                         tempQuantText = cell(length(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model})),1);
                         tempQuants = MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).(['Voxel_' num2str(gui.controls.act_x)])(gui.controls.Selected,:);
                         tempQuantText(idx_GABA) = table2cell(tempQuants(1,1));
                         idx_GABAp  = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{gui.quant.Selected.Model}),'GABAplus'));
                         tempQuantText(idx_GABAp) = table2cell(tempQuants(1,2));
                    end  
                    QuantText(2:end,q+1) = tempQuantText;
                else
                    QuantText(2:end,q+1) = table2cell(MRSCont.quantify.tables.(gui.quant.Names.Model{gui.quant.Selected.Model}).(gui.quant.Names.Quants{q}).(['Voxel_' num2str(gui.controls.act_x)])(gui.controls.Selected,:))';
                end
            end
            temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantText,...
                'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
             set( temp, 'BackgroundColor',gui.colormap.Background);
            delete(gui.Plot.quant.Children)
            set( temp, 'Parent', gui.Plot.quant ); %Update table
            set(gui.upperBox.quant.Info,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}]); %Update info Title
       else
           temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.quant.Names.Model{gui.quant.Selected.Model},[gui.controls.act_x gui.controls.act_y]); %Create figure
           ViewAxes = gca();
           
            delete(gui.Plot.quantHBox.Children(2).Children)
            set(ViewAxes.Children, 'Parent', gui.Plot.quantHBox.Children(2)); %Update plot
            set(gui.Plot.quantHBox.Children(2).Title, 'String', ViewAxes.Title.String) %Update title
            set(gui.Plot.quantHBox.Children(2),'Children',flipud(gui.Plot.quantHBox.Children(2).Children));
            set(gui.Plot.quantHBox.Children(2), 'XLim', ViewAxes.XLim) % Update Xlim
            set(gui.Plot.quantHBox.Children(2), 'YLim', ViewAxes.YLim) % Update Ylim
            set(gui.Plot.quantHBox.Children(2), 'Units', 'normalized');
            set(gui.Plot.quantHBox.Children(2), 'OuterPosition', [0.5,0.02,0.46,0.98])
            % Get rid of the Load figure
            close(temp);
            
            nominator = [];
            denominator = [];
            % Go through the checkbox list
            for idx = 1: gui.quant.Number.Model 
                for idx2 = 1 : length(gui.quant.Names.Metabs.(gui.quant.Names.Model{idx}))
                    if ~isempty(gui.Plot.MetabNodesNominator{idx,idx2}) 
                        if gui.Plot.MetabNodesNominator{idx,idx2}.Checked
                            nominator{end+1} =  gui.Plot.MetabNodesNominator{idx,idx2}.Name;
                        end
                    end
                    if ~isempty(gui.Plot.MetabNodesDenominator{idx,idx2})
                        if gui.Plot.MetabNodesDenominator{idx,idx2}.Checked
                            denominator{end+1} =  gui.Plot.MetabNodesDenominator{idx,idx2}.Name;
                        end
                    end
                end
            end
            if isempty(nominator)
                nominator{1} = 'NAA';
            end
            nominator_spec = [];
            denominator_spec = [];
            for idx = 1: length(gui.Plot.quantNodesNominator)
                if gui.Plot.quantNodesNominator{idx}.PartiallyChecked
                    nominator_spec{end+1} =  gui.Plot.quantNodesNominator{idx}.Name;
                end
                if gui.Plot.quantNodesDenominator{idx}.PartiallyChecked
                    denominator_spec{end+1} =  gui.Plot.quantNodesDenominator{idx}.Name;
                end
            end
            if ~isempty(denominator_spec)
                temp = osp_plotMRSImap(MRSCont, gui.controls.Selected, nominator_spec{1},denominator_spec{1}, nominator, denominator);
            else
                temp = osp_plotMRSImap(MRSCont, gui.controls.Selected, nominator_spec{1},[], nominator, denominator);
            end
            ViewAxes = gca();
            drawnow
            set( gui.Plot.quantMRSImap.Children, 'ColorData', ViewAxes.ColorData );
            set(gui.Plot.quantMRSImap.Children, 'Title', ViewAxes.Title{2}) %Update title
            set(gui.Plot.quantMRSImap.Children, 'Units', 'normalized');
            set(gui.Plot.quantMRSImap.Children, 'OuterPosition', [-0.15,0,1.2,0.98])
            close(temp)
            
            % If it is Multivoxel data we have to update the Voxel Position
            % window
            if MRSCont.flags.isMRSI 
                temp = osp_plotRawMRSIpos(MRSCont, 1, [gui.controls.act_y gui.controls.act_x]);
                ViewAxes = gca();
                drawnow
                set( gui.layout.LocPanel.Children,'ColorData', ViewAxes.ColorData );
                close(temp)
            end
           
        end
        setappdata(gui.figure,'MRSCont',MRSCont);  % Write MRSCont into hidden container in gui class
end