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
        gui.Plot.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(2);
        gui.Info.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2);
        gui.InfoText.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children;
        gui.Results.fit  = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1);
        gui.Results.FitTextAmpl = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1).Children(1).Children(1);
        gui.Results.FitTextNames = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1).Children(1).Children(2);
        
        gui.layout.EmptyFitPlot = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if gui.fit.Selected == 1 %Metabolite data?
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
        else if gui.fit.Selected == 2 %Reference data?
        StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            else %Water data?
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            end
        end
        set(gui.InfoText.fit, 'String',sprintf(StatText))
        % Update amplitudes for the fit results panel based on the files in the MRSCont (Raw Amplitudes or Water-scaled if ref or water supplied)
        if   ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w') %Is not concateneted or is reference/water fit 
            gui.fit.Style = gui.fit.Names{gui.fit.Selected};
            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
        else %Is concatenated and not water/reference
            gui.fit.Style = 'conc';
            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
        end
        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
        if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
            if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl) %Names and amplitudes
                    NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,gui.controls.Selected}.name{m} ': \n']];
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
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                end
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl)
                    NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,gui.controls.Selected}.name{m} ': \n']];
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
%%%3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'off' );
        temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,1,gui.fit.Names{gui.fit.Selected});
        ViewAxes = gca();
        delete(gui.Plot.fit.Children)
        set(ViewAxes.Children, 'Parent', gui.Plot.fit); %Update plot
        set(gui.Plot.fit.Title, 'String', ViewAxes.Title.String) %Update title
        set(gui.Plot.fit, 'XLim', ViewAxes.XLim) % Update Xlim
        set(gui.Plot.fit, 'YLim', ViewAxes.YLim) % Update Ylim
        % Get rid of the Load figure
        close(temp);
        set(gui.Info.fit,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}] );
        setappdata(gui.figure,'MRSCont',MRSCont); % Get MRSCont from hidden container in gui class
end