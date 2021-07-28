function osp_updateFitWindow(gui)
%% osp_updateFitWindow
%   This function updates the fit tab.
%
%
%   USAGE:
%       osp_updateFitWindow(gui);
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
        Selection = gui.fit.Names{gui.fit.Selected};
        gui.Plot.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(2);
        gui.upperBox.fit.Info = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2); 
        if (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x)
        end
        gui.controls.b_save_fitTab = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(1).Children;
        gui.InfoText.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children;
        gui.Results.fit  = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1);
        gui.Results.FitTextAmpl = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1).Children(1).Children(1);
        gui.Results.FitTextNames = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1).Children(1).Children(2);
        
        gui.layout.EmptyFitPlot = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w') %Is not concateneted or is reference/water fit 
            gui.fit.Style = gui.fit.Names{gui.fit.Selected};
        else %Is concatenated and not water/reference
            gui.fit.Style = 'conc';
        end
        if ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            ph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
            ph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
            if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                refShift = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                refFWHM = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM; 
                iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1; 
            end
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            ph0 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
            ph1 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
            if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                refShift = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                refFWHM = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM; 
                iniph0 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                iniph1 = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1; 
            end
        else
            RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            ph0 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
            ph1 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
            if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                refShift = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                refFWHM = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM; 
                iniph0 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                iniph1 = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1; 
            end
        end
      
        if  ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') %Metabolite data?
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f') 'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(MRSCont.fit.basisSet.nMets) '; Number of MM/lipids: ' num2str(MRSCont.fit.basisSet.nMM) ...
                        ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected}) '; initial ph0: ' num2str(iniph0,'%1.2f') 'deg; initial ph1: ' num2str(iniph1,'%1.2f') 'deg'];
        else if strcmp (Selection, 'ref') %Reference data?
        StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
            else %Water data?
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
            end
        end
        set(gui.upperBox.fit.Info.Children(2).Children, 'String',sprintf(StatText))
        % Update amplitudes for the fit results panel based on the files in the MRSCont (Raw Amplitudes or Water-scaled if ref or water supplied)
        if ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    if MRSCont.flags.hasRef
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl)
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(gui.Results.fit, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                end
            end
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet{1,gui.controls.act_x}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    if MRSCont.flags.hasRef
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,gui.controls.act_x}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl)
                        NameText = [NameText, [MRSCont.fit.resBasisSet{1,gui.controls.act_x}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(gui.Results.fit, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                end
            end    
        else
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    if MRSCont.flags.hasRef
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl)
                        try
                            NameText = [NameText, [MRSCont.fit.resBasisSet{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                        catch
                            NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,1}.name{m} ': \n']];
                        end
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(gui.Results.fit, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                end
            end 
        end
%%%3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'off' );
        if ~isfield(MRSCont.flags,'isPRIAM') && ~isfield(MRSCont.flags,'isMRSI') && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,1,Selection); %Create figure
        elseif  isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,gui.controls.act_x,Selection); %Create figure
        else
            temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y],Selection); %Create figure
        end

        ViewAxes = gca();
        delete(gui.Plot.fit.Children)
        set(ViewAxes.Children, 'Parent', gui.Plot.fit); %Update plot
        set(gui.Plot.fit.Title, 'String', ViewAxes.Title.String) %Update title
        set(gui.Plot.fit,'Children',flipud(gui.Plot.fit.Children));
        set(gui.Plot.fit, 'XLim', ViewAxes.XLim) % Update Xlim
        set(gui.Plot.fit, 'YLim', ViewAxes.YLim) % Update Ylim
        set(gui.Plot.fit, 'Units', 'normalized');
        set(gui.Plot.fit, 'OuterPosition', [0.17,0.02,0.75,0.98])
        % Get rid of the Load figure
        close(temp);
        
        % If it is Multivoxel data we have to update the Voxel Position
        % window
        if MRSCont.flags.isMRSI 
            temp = osp_plotRawMRSIpos(MRSCont, 1, [gui.controls.act_y gui.controls.act_x]);
            ViewAxes = gca();
            drawnow
            set( gui.layout.LocPanel.Children,'ColorData', ViewAxes.ColorData );
            close(temp)
        end
        set(gui.upperBox.fit.Info.Children(2),'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}] );
        set(gui.controls.b_save_fitTab,'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); % Get MRSCont from hidden container in gui class
end