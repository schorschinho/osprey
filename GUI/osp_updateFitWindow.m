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
        basis = gui.controls.act_z;
        subspectrum = gui.controls.act_y;
        experiment = gui.controls.act_x;

        
        switch MRSCont.opts.fit.method
            case 'LCModel'
                gui.Results.FitTextCRLB = gui.Results.fit{gui.fit.Selected}.Children.Children(1);
                gui.Results.FitTextAmpl = gui.Results.fit{gui.fit.Selected}.Children.Children(2);
                gui.Results.FitTextNames = gui.Results.fit{gui.fit.Selected}.Children.Children(3);
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(4),'String',gui.controls.act_z)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(5),'String',gui.controls.act_y)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(6),'String',gui.controls.act_x)
            case 'Osprey'
                gui.Results.FitTextAmpl = gui.Results.fit{gui.fit.Selected}.Children.Children(1);
                gui.Results.FitTextNames = gui.Results.fit{gui.fit.Selected}.Children.Children(2);
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(4),'String',gui.controls.act_z)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(5),'String',gui.controls.act_y)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(4).Children(1).Children.Children(6),'String',gui.controls.act_x)
            case 'Osprey_gLCM'
                gui.Results.FitTextCRLB = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(1).Children.Children(1);
                gui.Results.FitTextAmpl = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(1).Children.Children(2);
                gui.Results.FitTextNames = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(1).Children.Children(3);
                gui.controls.ModelStep = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(4);
                ModelStep = gui.controls.ModelStep.Value;
                % Selection = MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum}.Data.spec_name;
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(3).Children(4).Children(1).Children.Children(4),'String',gui.controls.act_z)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(3).Children(4).Children(1).Children.Children(5),'String',gui.controls.act_y)
                set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(3).Children(4).Children(1).Children.Children(6),'String',gui.controls.act_x)
        end

        gui.layout.EmptyFitPlot = 0;


%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w') %Is not concateneted or is reference/water fit
            gui.fit.Style = gui.fit.Names{gui.fit.Selected};
        else %Is concatenated and not water/reference
            gui.fit.Style = 'conc';
        end

        
        % For this visualization, we will have to make a few
        % distinctions upfront since the modeling algorithms (LCModel
        % vs. Osprey) do not always return the same kinds of data, or they
        % return them in different formats.

        switch MRSCont.opts.fit.method
            case 'LCModel'
                % Number of metabolites and lipid/MM basis functions
                basisNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.name;
                nLip    = sum(~cellfun(@isempty, strfind(basisNames, 'Lip')));
                nMM     = sum(~cellfun(@isempty, strfind(basisNames, 'MM')));
                nMMLip  = nLip + nMM;
                nMets   = length(basisNames) - nMMLip;
                nComb   = sum(~cellfun(@isempty, strfind(basisNames, '+')));
                % Where are the metabolite names stored?
                basisSetNames = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.name;
                % Smaller fonts for the results
                resultsFontSize = 9;
            case 'Osprey'
                % Where are the metabolite names stored?
                if strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')
                    basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.(gui.fit.Style){gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.(gui.fit.Style){gui.controls.Selected}.spectralwidth))]){1};
                    basisSetNames = basisSet.name;
                else if strcmp(gui.fit.Style, 'conc')
                        basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.spectralwidth))]){basis,1};
                        basisSetNames = basisSet.name;
                    else
                        basisSet = MRSCont.fit.resBasisSet.(gui.fit.Style).(['np_sw_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.sz(1))) '_' num2str(round(MRSCont.processed.metab{gui.controls.Selected}.spectralwidth))]){basis,1,subspectrum};
                        basisSetNames = basisSet.name;
                    end
                end
                % Number of metabolites and lipid/MM basis functions
                nMets   = basisSet.nMets;
                nMMLip  = basisSet.nMM;
                 % Larger fonts for the results
                resultsFontSize = 11;
            case 'Osprey_gLCM'
                % Number of metabolites and lipid/MM basis functions
                if (isfield(MRSCont.fit.results.(gui.fit.Style){basis,1,subspectrum,1}.Options{1,ModelStep},'paraIndirect'))              
                    basisSetNames =MRSCont.fit.results.(gui.fit.Style){basis,1,subspectrum,1}.BasisSets.names(MRSCont.fit.results.(gui.fit.Style){basis,1,subspectrum}.BasisSets.includeInFit(ModelStep,:)==1);
                    DisplayExperiment = experiment;
                    experiment = 1;
                else
                    basisSetNames =MRSCont.fit.results.(gui.fit.Style){basis,1,subspectrum,experiment}.BasisSets.names(MRSCont.fit.results.(gui.fit.Style){basis,1,subspectrum}.BasisSets.includeInFit(ModelStep,:)==1);
                    DisplayExperiment = 1;
                end
                nMMLip = sum(find(contains(basisSetNames,'MM') + contains(basisSetNames,'Lip')));
                nMets   = length(basisSetNames)-nMMLip;
                 % Larger fonts for the results
                resultsFontSize = 11;
                scale = MRSCont.fit.results.(gui.fit.Style){1,gui.controls.Selected}.scale;
        end

        if ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
             switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.h2oarea .* MRSCont.fit.scale{1,gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{1,gui.controls.Selected};
                        CRLB    = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected,subspectrum}.CRLB;
                    end
                case 'Osprey'
                    RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{1,gui.controls.Selected};
                case 'Osprey_gLCM'
                    RawAmpl = MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.Model{ModelStep}.parsOut.metAmpl(DisplayExperiment,:) .* MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.scale;
                    CRLB    = MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.Model{ModelStep}.CRLB{1,:};

             end
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.h2oarea .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.CRLB;
                    end
                case 'Osprey_gLCM'
                    RawAmpl = MRSCont.fit.results{1,gui.controls.act_x}.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum,experiment}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                    CRLB    = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.Model{end}.CRLB{1,:};
            end
        else
            switch MRSCont.opts.fit.method
                case 'LCModel'
                    if strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.h2oarea .* MRSCont.fit.scale{gui.controls.Selected};
                    else
                        RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                        CRLB    = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{gui.controls.Selected,subspectrum}.CRLB;
                    end
                case 'Osprey'
                    RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style).fitParams{basis,gui.controls.Selected,subspectrum}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
                case 'Osprey_gLCM'
                    RawAmpl = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style){gui.controls.Selected,end,subspectrum,experiment}.Model{ModelStep}.parsOut.metAmpl .* MRSCont.fit.results.(gui.fit.Style){gui.controls.Selected,end}.scale;
                    CRLB    = MRSCont.fit.results{gui.controls.act_x,gui.controls.act_y}.(gui.fit.Style){gui.controls.Selected,end,subspectrum,experiment}.Model{ModelStep}.CRLB{1,:};
            end
        end

        if  ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') %Metabolite data?
            StatText = [ 'Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Selected subspecs: ' Selection ];
        else if strcmp (Selection, 'ref') %Reference data?
                StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Selected subspecs: ' Selection];
            else %Is water data
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Selected subspecs: ' Selection];
            end
        end
        set(gui.upperBox.fit.Info{gui.fit.Selected}.Children, 'String',sprintf(StatText))
        % Update amplitudes for the fit results panel based on the files in the MRSCont (Raw Amplitudes or Water-scaled if ref or water supplied)
        if ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                NameText    = [''];
                RawAmplText = [''];
                CRLBText    = [''];
                for m = 1 : length(RawAmpl) %Names and Amplitudes
                    NameText = [NameText, [basisSetNames{m} ' \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    if strcmp(MRSCont.opts.fit.method, 'LCModel') || strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                        if isinf(CRLB(m))
                            CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '\n']]; 
                        else if CRLB(m) > 999
                                CRLBText = [CRLBText, [num2str(Inf, '%1.3g') '\n']]; 
                            else
                                CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '%%\n']];
                            end
                        end
                        
                    end

                end
            else %Water/reference fit but this should never happen in this loop
                NameText = ['Water: ' ];
                RawAmplText = [num2str(RawAmpl,'%1.2e')];
            end
                set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
                set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                if strcmp(MRSCont.opts.fit.method, 'LCModel') || strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                    set(gui.Results.FitTextCRLB, 'String',sprintf(CRLBText));
                end
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit?
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                        end
                        case 'LCModel'
                        case 'Osprey_gLCM'
                            if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                                RawAmpl = RawAmpl ./ (sum(MRSCont.fit.results.ref{1,gui.controls.Selected}.Model{1, 1}.parsOut.metAmpl(DisplayExperiment,:)) .* MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.scale);
                            else
                                RawAmpl = RawAmpl ./ (sum(MRSCont.fit.results.w{1,gui.controls.Selected}.Model{1, 1}.parsOut.metAmpl(DisplayExperiment,:)) .* MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.scale);
                            end
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    CRLBText    = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [basisSetNames{m} ' \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        if strcmp(MRSCont.opts.fit.method, 'LCModel') || strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                            if isinf(CRLB(m))
                                CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '\n']]; 
                            else if CRLB(m) > 999
                                    CRLBText = [CRLBText, [num2str(Inf, '%1.3g') '\n']]; 
                                else
                                    CRLBText = [CRLBText, [num2str(round(CRLB(m),1), '%1.3g') '%%\n']];
                                end
                            end
                        end

                    end
                    set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                    if strcmp(MRSCont.opts.fit.method, 'LCModel') || strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                        set(gui.Results.FitTextCRLB, 'String',sprintf(CRLBText));
                    end
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
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
                set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
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
                    set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
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
                set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
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
                    set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Water Ratio']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                else %Water fit
                   NameText = ['Water: \t'];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit{gui.fit.Selected}, 'Title', ['Raw Amplitudes']);
                    set(gui.Results.FitTextNames, 'String',sprintf(NameText));
                    set(gui.Results.FitTextAmpl, 'String',sprintf(RawAmplText));
                end
            end
        end
%%%3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'on' );
        if ~strcmp(MRSCont.opts.fit.method,'Osprey_gLCM')
            if  ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y gui.controls.act_z],Selection); %Create figure
            elseif  isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,gui.controls.act_x,Selection); %Create figure
            else
                temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,[gui.controls.act_x gui.controls.act_y],Selection); %Create figure
            end
        else
            switch gui.overview.Selected.ModelPlot
                case 1
                    MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.plotFit1D(0,ModelStep,DisplayExperiment);
                case 2
                    MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.plotFit1DStack(0,ModelStep,DisplayExperiment);
                case 3
                    if size(MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.Data.fids) > 1
                        MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.plotFit3D(0,ModelStep);
                    else
                        MRSCont.fit.results.(gui.fit.Style){basis,gui.controls.Selected,subspectrum,experiment}.plotFit1D(0,ModelStep,DisplayExperiment);
                    end
            end         
        end

        ViewAxes = gca();
        tempXLim = ViewAxes.XLim;
        tempYLim = ViewAxes.YLim;
        tempZLim = ViewAxes.ZLim;
        delete(gui.Plot.fit{gui.fit.Selected}.Children(2).Children)
        set(ViewAxes.Children, 'Parent', gui.Plot.fit{gui.fit.Selected}.Children(2)); %Update plot
        set(gui.Plot.fit{gui.fit.Selected}.Children(2).Title, 'String', ViewAxes.Title.String) %Update title
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'View', ViewAxes.View) % View     
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'XLim', tempXLim) % Update Xlim
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'YLim', tempYLim) % Update Ylim
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'ZLim', tempZLim) % Update Zlim
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'XTick', ViewAxes.XTick) % Update XTick
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'YTick', ViewAxes.YTick) % Update YTick
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'ZTick', ViewAxes.ZTick) % Update ZTick
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'XTickLabel', ViewAxes.XTickLabel) % Update XTickLabel
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'YTickLabel', ViewAxes.YTickLabel) % Update YTickLabel
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'ZTickLabel', ViewAxes.ZTickLabel) % Update ZTickLabel
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'XColor', ViewAxes.XColor) % Update XColor
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'YColor', gui.colormap.Background) % Update YColor
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'ZColor', gui.colormap.Background) % Update ZColor   
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'XDir', ViewAxes.XDir) % XDir
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'YDir', ViewAxes.YDir) % YDir
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'ZDir', ViewAxes.ZDir) % ZDir
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'Units', 'normalized');
        set(gui.Plot.fit{gui.fit.Selected}.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
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
        set(gui.upperBox.fit.Info{gui.fit.Selected},'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}] );
        set(gui.controls.b_save_fitTab{gui.fit.Selected},'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); % Get MRSCont from hidden container in gui class
end
