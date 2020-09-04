classdef CreateOspreyJob_app < matlab.apps.AppBase
%%  CreateOspreyJob_app
%   This app creates a windwo to interactivley create an OspreyJob.m file.
%   You can specify fit options, spectra, anatomical images, and included
%   metabolites during this process. It is called from the Osprey GUI
%   window.
%
%   AUTHOR:
%       Dr. Peter Van Schuerbeek (UZ Brussel (VUB), 2020-09-03)
%
%   CREDITS:
%       The inital code is was created and kindly shared by Dr. Peter Van
%       Schuerbeek. We adapted this version based on his contribution.
%
%   HISTORY:
%       2020-09-03: First version of the code.
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SpecifySequenceInformationPanel  matlab.ui.container.Panel
        SequenceTypeDropDownLabel       matlab.ui.control.Label
        SequenceTypeDropDown            matlab.ui.control.DropDown
        EditingTargetsDropDownLabel     matlab.ui.control.Label
        EditingTargetsDropDown          matlab.ui.control.DropDown
        TextArea                        matlab.ui.control.TextArea
        SpecifyDataHandlingandModelingOptionsPanel  matlab.ui.container.Panel
        SaveLCMCheckBox                 matlab.ui.control.CheckBox
        SaveJMRUICheckBox               matlab.ui.control.CheckBox
        SaveVendorCheckBox              matlab.ui.control.CheckBox
        FittingAlgorithmDropDownLabel   matlab.ui.control.Label
        FittingAlgorithmDropDown        matlab.ui.control.DropDown
        IncludedMetabolitesDropDownLabel  matlab.ui.control.Label
        IncludedMetabolitesDropDown     matlab.ui.control.DropDown
        FittingStyleDropDownLabel       matlab.ui.control.Label
        FittingStyleDropDown            matlab.ui.control.DropDown
        MRSFitRangeppmEditFieldLabel    matlab.ui.control.Label
        MRSFitRangeppmEditField         matlab.ui.control.EditField
        WaterFitRangeppmEditFieldLabel  matlab.ui.control.Label
        WaterFitRangeppmEditField       matlab.ui.control.EditField
        BaselineknotspacingppmEditFieldLabel  matlab.ui.control.Label
        BaselineknotspacingppmEditField  matlab.ui.control.NumericEditField
        AddMMandLipbasisfunctionstofitCheckBox  matlab.ui.control.CheckBox
        SelectedMetabolitesPanel        matlab.ui.container.Panel
        AlaCheckBox                     matlab.ui.control.CheckBox
        AscCheckBox                     matlab.ui.control.CheckBox
        AspCheckBox                     matlab.ui.control.CheckBox
        bHBCheckBox                     matlab.ui.control.CheckBox
        bHGCheckBox                     matlab.ui.control.CheckBox
        CitCheckBox                     matlab.ui.control.CheckBox
        CrCheckBox                      matlab.ui.control.CheckBox
        CrCH2CheckBox                   matlab.ui.control.CheckBox
        EtOHCheckBox                    matlab.ui.control.CheckBox
        GABACheckBox                    matlab.ui.control.CheckBox
        GPCCheckBox                     matlab.ui.control.CheckBox
        GSHCheckBox                     matlab.ui.control.CheckBox
        GlcCheckBox                     matlab.ui.control.CheckBox
        GlnCheckBox                     matlab.ui.control.CheckBox
        GluCheckBox                     matlab.ui.control.CheckBox
        GlyCheckBox                     matlab.ui.control.CheckBox
        H2OCheckBox                     matlab.ui.control.CheckBox
        InsCheckBox                     matlab.ui.control.CheckBox
        LacCheckBox                     matlab.ui.control.CheckBox
        NAACheckBox                     matlab.ui.control.CheckBox
        NAAGCheckBox                    matlab.ui.control.CheckBox
        PChCheckBox                     matlab.ui.control.CheckBox
        PCrCheckBox                     matlab.ui.control.CheckBox
        PECheckBox                      matlab.ui.control.CheckBox
        PhenylCheckBox                  matlab.ui.control.CheckBox
        ScylloCheckBox                  matlab.ui.control.CheckBox
        SerCheckBox                     matlab.ui.control.CheckBox
        TauCheckBox                     matlab.ui.control.CheckBox
        TyrosCheckBox                   matlab.ui.control.CheckBox
        MM09CheckBox                    matlab.ui.control.CheckBox
        MM12CheckBox                    matlab.ui.control.CheckBox
        MM14CheckBox                    matlab.ui.control.CheckBox
        MM17CheckBox                    matlab.ui.control.CheckBox
        MM20CheckBox                    matlab.ui.control.CheckBox
        Lip09CheckBox                   matlab.ui.control.CheckBox
        Lip13CheckBox                   matlab.ui.control.CheckBox
        Lip20CheckBox                   matlab.ui.control.CheckBox
        SpecifyMRSandStructuralImagingFilesPanel  matlab.ui.container.Panel
        MRSDataButton                   matlab.ui.control.Button
        H2OReferenceButton              matlab.ui.control.Button
        H2OShortTEButton                matlab.ui.control.Button
        MetaboliteNulledButton          matlab.ui.control.Button
        T1DataniftiniiButton            matlab.ui.control.Button
        NumberofdatasetsEditFieldLabel  matlab.ui.control.Label
        NumberofdatasetsEditField       matlab.ui.control.NumericEditField
        MRSDataText                     matlab.ui.control.TextArea
        H2OShortTEText                  matlab.ui.control.TextArea
        H2OReferenceText                matlab.ui.control.TextArea
        MetaboliteNulledText            matlab.ui.control.TextArea
        T1DataText                      matlab.ui.control.TextArea
        SpecifyOutputFolderPanel        matlab.ui.container.Panel
        OutputFolderEditField           matlab.ui.control.EditField
        OutputFolderButton              matlab.ui.control.Button
        JobNameEditFieldLabel           matlab.ui.control.Label
        JobNameEditField                matlab.ui.control.EditField
        CANCELButton                    matlab.ui.control.Button
        CREATEJOBButton                 matlab.ui.control.Button
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: SequenceTypeDropDown
        function SequenceTypeDropDownValueChanged(app, event)
            value = app.SequenceTypeDropDown.Value;
            
            switch value
                case 'unedited'
                    app.EditingTargetsDropDown.Items = {'none'};
                    app.EditingTargetsDropDown.Value = 'none';
                    
                    app.FittingStyleDropDown.Enable = 'Off';
                    app.FittingStyleDropDownLabel.Enable = 'Off';
                case 'MEGA'
                    app.EditingTargetsDropDown.Items = {'GABA','GSH','Lac'};
                    app.EditingTargetsDropDown.Value = 'GABA';
                    
                    app.FittingStyleDropDown.Enable = 'On';
                    app.FittingStyleDropDownLabel.Enable = 'On';
                case 'HERMES'
                    app.EditingTargetsDropDown.Items = {'GABA, GSH','GABA, GSH, EtOH'};
                    app.EditingTargetsDropDown.Value = 'GABA, GSH';
                    
                    app.FittingStyleDropDown.Enable = 'On';
                    app.FittingStyleDropDownLabel.Enable = 'On';
                case 'HERCULES'
                    app.EditingTargetsDropDown.Items = {'HERCULES1','HERCULES2','HERCULES3'};
                    app.EditingTargetsDropDown.Value = 'HERCULES1';
                    
                    app.FittingStyleDropDown.Enable = 'On';
                    app.FittingStyleDropDownLabel.Enable = 'On';
            end
        end

        % Value changed function: IncludedMetabolitesDropDown
        function IncludedMetabolitesDropDownValueChanged(app, event)
            value = app.IncludedMetabolitesDropDown.Value;
            
            switch value
                case 'default'
                    app.AlaCheckBox.Value = false;
                    app.AlaCheckBox.Enable = 'Off';
                    app.AscCheckBox.Value = true;
                    app.AscCheckBox.Enable = 'Off';
                    app.AspCheckBox.Value = true;
                    app.AspCheckBox.Enable = 'Off';
                    app.bHBCheckBox.Value = false;
                    app.bHBCheckBox.Enable = 'Off';
                    app.bHGCheckBox.Value = false;
                    app.bHGCheckBox.Enable = 'Off';
                    app.CitCheckBox.Value = false;
                    app.CitCheckBox.Enable = 'Off';
                    app.CrCheckBox.Value = true;
                    app.CrCheckBox.Enable = 'Off';
                    app.CrCH2CheckBox.Value = true;
                    app.CrCH2CheckBox.Enable = 'Off';
                    app.EtOHCheckBox.Value = false;
                    app.EtOHCheckBox.Enable = 'Off';
                    app.GABACheckBox.Value = true;
                    app.GABACheckBox.Enable = 'Off';
                    app.GPCCheckBox.Value = true;
                    app.GPCCheckBox.Enable = 'Off';
                    app.GSHCheckBox.Value = true;
                    app.GSHCheckBox.Enable = 'Off';
                    app.GlcCheckBox.Value = false;
                    app.GlcCheckBox.Enable = 'Off';
                    app.GlnCheckBox.Value = true;
                    app.GlnCheckBox.Enable = 'Off';
                    app.GluCheckBox.Value = true;
                    app.GluCheckBox.Enable = 'Off';
                    app.GlyCheckBox.Value = false;
                    app.GlyCheckBox.Enable = 'Off';
                    app.H2OCheckBox.Value = true;
                    app.H2OCheckBox.Enable = 'Off';
                    app.InsCheckBox.Value = true;
                    app.InsCheckBox.Enable = 'Off';
                    app.LacCheckBox.Value = true;
                    app.LacCheckBox.Enable = 'Off';
                    app.NAACheckBox.Value = true;
                    app.NAACheckBox.Enable = 'Off';
                    app.NAAGCheckBox.Value = true;
                    app.NAAGCheckBox.Enable = 'Off';
                    app.PChCheckBox.Value = true;
                    app.PChCheckBox.Enable = 'Off';
                    app.PCrCheckBox.Value = true;
                    app.PCrCheckBox.Enable = 'Off';
                    app.PECheckBox.Value = true;
                    app.PECheckBox.Enable = 'Off';
                    app.PhenylCheckBox.Value = false;
                    app.PhenylCheckBox.Enable = 'Off';
                    app.ScylloCheckBox.Value = true;
                    app.ScylloCheckBox.Enable = 'Off';
                    app.SerCheckBox.Value = false;
                    app.SerCheckBox.Enable = 'Off';
                    app.TauCheckBox.Value = true;
                    app.TauCheckBox.Enable = 'Off';
                    app.TyrosCheckBox.Value = false;
                    app.TyrosCheckBox.Enable = 'Off';
                    app.MM09CheckBox.Value = true;
                    app.MM09CheckBox.Enable = 'Off';
                    app.MM12CheckBox.Value = true;
                    app.MM12CheckBox.Enable = 'Off';
                    app.MM14CheckBox.Value = true;
                    app.MM14CheckBox.Enable = 'Off';
                    app.MM17CheckBox.Value = true;
                    app.MM17CheckBox.Enable = 'Off';
                    app.MM20CheckBox.Value = true;
                    app.MM20CheckBox.Enable = 'Off';
                    app.Lip09CheckBox.Value = true;
                    app.Lip09CheckBox.Enable = 'Off';
                    app.Lip13CheckBox.Value = true;
                    app.Lip13CheckBox.Enable = 'Off';
                    app.Lip20CheckBox.Value = true;
                    app.Lip20CheckBox.Enable = 'Off';
                case 'full'
                    app.AlaCheckBox.Value = true;
                    app.AlaCheckBox.Enable = 'Off';
                    app.AscCheckBox.Value = true;
                    app.AscCheckBox.Enable = 'Off';
                    app.AspCheckBox.Value = true;
                    app.AspCheckBox.Enable = 'Off';
                    app.bHBCheckBox.Value = true;
                    app.bHBCheckBox.Enable = 'Off';
                    app.bHGCheckBox.Value = true;
                    app.bHGCheckBox.Enable = 'Off';
                    app.CitCheckBox.Value = true;
                    app.CitCheckBox.Enable = 'Off';
                    app.CrCheckBox.Value = true;
                    app.CrCheckBox.Enable = 'Off';
                    app.CrCH2CheckBox.Value = true;
                    app.CrCH2CheckBox.Enable = 'Off';
                    app.EtOHCheckBox.Value = true;
                    app.EtOHCheckBox.Enable = 'Off';
                    app.GABACheckBox.Value = true;
                    app.GABACheckBox.Enable = 'Off';
                    app.GPCCheckBox.Value = true;
                    app.GPCCheckBox.Enable = 'Off';
                    app.GSHCheckBox.Value = true;
                    app.GSHCheckBox.Enable = 'Off';
                    app.GlcCheckBox.Value = true;
                    app.GlcCheckBox.Enable = 'Off';
                    app.GlnCheckBox.Value = true;
                    app.GlnCheckBox.Enable = 'Off';
                    app.GluCheckBox.Value = true;
                    app.GluCheckBox.Enable = 'Off';
                    app.GlyCheckBox.Value = true;
                    app.GlyCheckBox.Enable = 'Off';
                    app.H2OCheckBox.Value = true;
                    app.H2OCheckBox.Enable = 'Off';
                    app.InsCheckBox.Value = true;
                    app.InsCheckBox.Enable = 'Off';
                    app.LacCheckBox.Value = true;
                    app.LacCheckBox.Enable = 'Off';
                    app.NAACheckBox.Value = true;
                    app.NAACheckBox.Enable = 'Off';
                    app.NAAGCheckBox.Value = true;
                    app.NAAGCheckBox.Enable = 'Off';
                    app.PChCheckBox.Value = true;
                    app.PChCheckBox.Enable = 'Off';
                    app.PCrCheckBox.Value = true;
                    app.PCrCheckBox.Enable = 'Off';
                    app.PECheckBox.Value = true;
                    app.PECheckBox.Enable = 'Off';
                    app.PhenylCheckBox.Value = true;
                    app.PhenylCheckBox.Enable = 'Off';
                    app.ScylloCheckBox.Value = true;
                    app.ScylloCheckBox.Enable = 'Off';
                    app.SerCheckBox.Value = true;
                    app.SerCheckBox.Enable = 'Off';
                    app.TauCheckBox.Value = true;
                    app.TauCheckBox.Enable = 'Off';
                    app.TyrosCheckBox.Value = true;
                    app.TyrosCheckBox.Enable = 'Off';
                    app.MM09CheckBox.Value = true;
                    app.MM09CheckBox.Enable = 'Off';
                    app.MM12CheckBox.Value = true;
                    app.MM12CheckBox.Enable = 'Off';
                    app.MM14CheckBox.Value = true;
                    app.MM14CheckBox.Enable = 'Off';
                    app.MM17CheckBox.Value = true;
                    app.MM17CheckBox.Enable = 'Off';
                    app.MM20CheckBox.Value = true;
                    app.MM20CheckBox.Enable = 'Off';
                    app.Lip09CheckBox.Value = true;
                    app.Lip09CheckBox.Enable = 'Off';
                    app.Lip13CheckBox.Value = true;
                    app.Lip13CheckBox.Enable = 'Off';
                    app.Lip20CheckBox.Value = true;
                    app.Lip20CheckBox.Enable = 'Off';
                case 'custom'
                    app.AlaCheckBox.Enable = 'On';
                    app.AscCheckBox.Enable = 'On';
                    app.AspCheckBox.Enable = 'On';
                    app.bHBCheckBox.Enable = 'On';
                    app.bHGCheckBox.Enable = 'On';
                    app.CitCheckBox.Enable = 'On';
                    app.CrCheckBox.Enable = 'Off';
                    app.CrCH2CheckBox.Enable = 'On';
                    app.EtOHCheckBox.Enable = 'On';
                    app.GABACheckBox.Enable = 'On';
                    app.GPCCheckBox.Enable = 'Off';
                    app.GSHCheckBox.Enable = 'On';
                    app.GlcCheckBox.Enable = 'On';
                    app.GlnCheckBox.Enable = 'On';
                    app.GluCheckBox.Enable = 'Off';
                    app.GlyCheckBox.Enable = 'On';
                    app.H2OCheckBox.Enable = 'On';
                    app.InsCheckBox.Enable = 'Off';
                    app.LacCheckBox.Enable = 'On';
                    app.NAACheckBox.Enable = 'Off';
                    app.NAAGCheckBox.Enable = 'On';
                    app.PChCheckBox.Enable = 'On';
                    app.PCrCheckBox.Enable = 'On';
                    app.PECheckBox.Enable = 'On';
                    app.PhenylCheckBox.Enable = 'On';
                    app.ScylloCheckBox.Enable = 'On';
                    app.SerCheckBox.Enable = 'On';
                    app.TauCheckBox.Enable = 'On';
                    app.TyrosCheckBox.Enable = 'On';
                    app.MM09CheckBox.Enable = 'On';
                    app.MM12CheckBox.Enable = 'On';
                    app.MM14CheckBox.Enable = 'On';
                    app.MM17CheckBox.Enable = 'On';
                    app.MM20CheckBox.Enable = 'On';
                    app.Lip09CheckBox.Enable = 'On';
                    app.Lip13CheckBox.Enable = 'On';
                    app.Lip20CheckBox.Enable = 'On';
            end
        end

        % Button pushed function: MRSDataButton
        function MRSDataButtonPushed(app, event)
            info = 'Please select the metabolite spectra to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            mrsfiles = spm_select(ndata,'any',info,{},pwd,'.*','1');

            [~,file_basename,file_exten]=fileparts(mrsfiles(1,:));
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} mrsfiles(i,:)};
            end
            
            app.MRSDataText.Value = filelist;
            
            if strcmp(file_exten,'.7')
                app.H2OReferenceButton.Enable = 'Off';
            else
                app.H2OReferenceButton.Enable = 'On';
            end
            
            app.H2OShortTEButton.Enable = 'On';
            app.MetaboliteNulledButton.Enable = 'On';
            app.T1DataniftiniiButton.Enable = 'On';
            app.NumberofdatasetsEditField.Enable = 'Off';
            movegui(app.UIFigure,'onscreen')
        end

        % Button pushed function: H2OReferenceButton
        function H2OReferenceButtonPushed(app, event)
            info = 'Please select the water reference file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            dir = app.MRSDataText.Value{1};
            [path,~,~]=fileparts(dir);
            
            SepFileList =  split(path, filesep);
            npath = [];
            for s = 1 : length(SepFileList)-1
                npath = [npath SepFileList{s} filesep];
            end
            
            try
                h2oreffiles = spm_select(ndata,'any',info,{},npath,'.*','1');
            catch
                h2oreffiles = spm_select(ndata,'any',info,{},path,'.*','1');
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2oreffiles(i,:)};
            end
            
            app.H2OReferenceText.Value = filelist;
            movegui(app.UIFigure,'onscreen')
        end

        % Button pushed function: H2OShortTEButton
        function H2OShortTEButtonPushed(app, event)
            info = 'Please select the water short-TE file to read';
            [fname,pathname]=uigetfile('*.*',info);
            
            ndata = app.NumberofdatasetsEditField.Value;
            dir = app.MRSDataText.Value{1};
            [path,~,~]=fileparts(dir);
            SepFileList =  split(path, filesep);
            npath = [];
            for s = 1 : length(SepFileList)-1
                npath = [npath SepFileList{s} filesep];
            end
            
            try           
                h2ostefiles = spm_select(ndata,'any',info,{},npath,'.*','1');
            catch
                h2ostefiles = spm_select(ndata,'any',info,{},path,'.*','1');
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2ostefiles(i,:)};
            end
            
            app.H2OShortTEText.Value = filelist;
            movegui(app.UIFigure,'onscreen')
        end

        % Button pushed function: MetaboliteNulledButton
        function MetaboliteNulledButtonPushed(app, event)
            info = 'Please select the metabolite-nulled file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            dir = app.MRSDataText.Value{1};
            [path,~,~]=fileparts(dir);
            SepFileList =  split(path, filesep);
            npath = [];
            for s = 1 : length(SepFileList)-1
                npath = [npath SepFileList{s} filesep];
            end
            
            try 
                metnulfiles = spm_select(ndata,'any',info,{},npath,'.*','1');
            catch
                metnulfiles = spm_select(ndata,'any',info,{},path,'.*','1');
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} metnulfiles(i,:)};
            end
            
            app.MetaboliteNulledText.Value = filelist;
            movegui(app.UIFigure,'onscreen')
        end

        % Button pushed function: T1DataniftiniiButton
        function T1DataniftiniiButtonPushed(app, event)
            info = 'Please select the T1 anatomical file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            dir = app.MRSDataText.Value{1};
            [path,~,~]=fileparts(dir);
            SepFileList =  split(path, filesep);
            npath = [];
            for s = 1 : length(SepFileList)-2
                npath = [npath SepFileList{s} filesep];
            end
            
            try            
                t1imfiles = spm_select(ndata,'any',info,{},npath,'.*','1');
            catch
                t1imfiles = spm_select(ndata,'any',info,{},path,'.*','1');
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} t1imfiles(i,:)};
            end
            
            app.T1DataText.Value = filelist;
            movegui(app.UIFigure,'onscreen')
        end

        % Button pushed function: OutputFolderButton
        function OutputFolderButtonPushed(app, event)
            info = 'Please select the output folder';
            pathname=uigetdir('*.*',info);
            
            app.OutputFolderEditField.Value = pathname;
        end

        % Button pushed function: CANCELButton
        function CANCELButtonPushed(app, event)
            delete(app.UIFigure)
        end

        % Button pushed function: CREATEJOBButton
        function CREATEJOBButtonPushed(app, event)
            jobm = osp_create_job_file(app);
            
            delete(app.UIFigure)
            
            MRSCont = OspreyJob(jobm,1);
            idx = 1; %gui.colors.Value;

            %Default colormap
            colormap.Background = [255/255 254/255 254/255];
            colormap.LightAccent = [110/255 136/255 164/255];
            colormap.Foreground = [11/255 71/255 111/255];
            colormap.Accent = [254/255 186/255 47/255];

            MRSCont.colormap = colormap;
            MRSCont.colormapidx = idx;
            MRSCont.flags.isToolChecked = 1;
            OspreyGUI(MRSCont);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 715 600];
            app.UIFigure.Name = 'Interactive OspreyJob.m Creator';
            app.UIFigure.Color = [255/255 254/255 254/255];

            % Create SpecifySequenceInformationPanel
            app.SpecifySequenceInformationPanel = uipanel(app.UIFigure);
            app.SpecifySequenceInformationPanel.ForegroundColor = [11/255 71/255 111/255];
            app.SpecifySequenceInformationPanel.BackgroundColor = [255/255 254/255 254/255];
            app.SpecifySequenceInformationPanel.Title = '1. Specify Sequence Information';
            app.SpecifySequenceInformationPanel.FontWeight = 'bold';
            app.SpecifySequenceInformationPanel.FontSize = 15;
            app.SpecifySequenceInformationPanel.FontName = 'Arial';
            app.SpecifySequenceInformationPanel.Position = [5 495 350 100];

            % Create SequenceTypeDropDownLabel
            app.SequenceTypeDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.SequenceTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SequenceTypeDropDownLabel.Position = [28 45 88 22];
            app.SequenceTypeDropDownLabel.Text = 'Sequence Type';
            app.SequenceTypeDropDownLabel.FontColor = [11/255 71/255 111/255];

            % Create SequenceTypeDropDown
            app.SequenceTypeDropDown = uidropdown(app.SpecifySequenceInformationPanel,'FontColor',[11/255 71/255 111/255]);
            app.SequenceTypeDropDown.BackgroundColor = [255/255 254/255 254/255];
            app.SequenceTypeDropDown.Items = {'unedited', 'MEGA', 'HERMES', 'HERCULES'};
            app.SequenceTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SequenceTypeDropDownValueChanged, true);
            app.SequenceTypeDropDown.Position = [131 45 129 22];
            app.SequenceTypeDropDown.Value = 'unedited';

            % Create EditingTargetsDropDownLabel
            app.EditingTargetsDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.EditingTargetsDropDownLabel.HorizontalAlignment = 'right';
            app.EditingTargetsDropDownLabel.Position = [28 12 85 22];
            app.EditingTargetsDropDownLabel.Text = 'Editing Targets';
            app.EditingTargetsDropDownLabel.FontColor = [11/255 71/255 111/255];

            % Create EditingTargetsDropDown
            app.EditingTargetsDropDown = uidropdown(app.SpecifySequenceInformationPanel,'FontColor',[11/255 71/255 111/255]);
            app.EditingTargetsDropDown.BackgroundColor = [255/255 254/255 254/255];
            app.EditingTargetsDropDown.Items = {'none'};
            app.EditingTargetsDropDown.Position = [131 12 129 22];
            app.EditingTargetsDropDown.Value = 'none';

            % Create SpecifyDataHandlingandModelingOptionsPanel
            app.SpecifyDataHandlingandModelingOptionsPanel = uipanel(app.UIFigure);
            app.SpecifyDataHandlingandModelingOptionsPanel.ForegroundColor = [11/255 71/255 115/255];
            app.SpecifyDataHandlingandModelingOptionsPanel.BackgroundColor = [255/255 254/255 254/255];
            app.SpecifyDataHandlingandModelingOptionsPanel.Title = '2. Specify Data Handling and Modeling Options';
            app.SpecifyDataHandlingandModelingOptionsPanel.FontWeight = 'bold';
            app.SpecifyDataHandlingandModelingOptionsPanel.FontSize = 15;
            app.SpecifyDataHandlingandModelingOptionsPanel.Position = [5 170 350 313];

            % Create SaveLCMCheckBox
            app.SaveLCMCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveLCMCheckBox.Text = 'Save LCM';
            app.SaveLCMCheckBox.Position = [21 253 78 22];
            app.SaveLCMCheckBox.FontColor = [11/255 71/255 111/255];
            

            % Create SaveJMRUICheckBox
            app.SaveJMRUICheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveJMRUICheckBox.Text = 'Save JMRUI';
            app.SaveJMRUICheckBox.Position = [124 253 89 22];
            app.SaveJMRUICheckBox.FontColor = [11/255 71/255 111/255];

            % Create SaveVendorCheckBox
            app.SaveVendorCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveVendorCheckBox.Text = 'Save Vendor';
            app.SaveVendorCheckBox.Position = [238 253 90 22];
            app.SaveVendorCheckBox.FontColor = [11/255 71/255 111/255];

            % Create FittingAlgorithmDropDownLabel
            app.FittingAlgorithmDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingAlgorithmDropDownLabel.HorizontalAlignment = 'right';
            app.FittingAlgorithmDropDownLabel.Position = [28 212 94 22];
            app.FittingAlgorithmDropDownLabel.Text = 'Fitting Algorithm';
            app.FittingAlgorithmDropDownLabel.FontColor = [11/255 71/255 111/255];

            % Create FittingAlgorithmDropDown
            app.FittingAlgorithmDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel,'FontColor',[11/255 71/255 111/255]);
            app.FittingAlgorithmDropDown.BackgroundColor = [255/255 254/255 254/255];
            app.FittingAlgorithmDropDown.Items = {'Osprey'};
            app.FittingAlgorithmDropDown.Position = [161 212 129 22];
            app.FittingAlgorithmDropDown.Value = 'Osprey';

            % Create IncludedMetabolitesDropDownLabel
            app.IncludedMetabolitesDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.IncludedMetabolitesDropDownLabel.HorizontalAlignment = 'right';
            app.IncludedMetabolitesDropDownLabel.Position = [28 179 118 22];
            app.IncludedMetabolitesDropDownLabel.Text = 'Included Metabolites';
            app.IncludedMetabolitesDropDownLabel.FontColor = [11/255 71/255 111/255];

            % Create IncludedMetabolitesDropDown
            app.IncludedMetabolitesDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel,'FontColor',[11/255 71/255 111/255]);
            app.IncludedMetabolitesDropDown.BackgroundColor = [255/255 254/255 254/255];
            app.IncludedMetabolitesDropDown.Items = {'default', 'full', 'custom'};
            app.IncludedMetabolitesDropDown.ValueChangedFcn = createCallbackFcn(app, @IncludedMetabolitesDropDownValueChanged, true);
            app.IncludedMetabolitesDropDown.Position = [161 179 129 22];
            app.IncludedMetabolitesDropDown.Value = 'default';

            % Create FittingStyleDropDownLabel
            app.FittingStyleDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingStyleDropDownLabel.HorizontalAlignment = 'right';
            app.FittingStyleDropDownLabel.Enable = 'off';
            app.FittingStyleDropDownLabel.Position = [28 145 69 22];
            app.FittingStyleDropDownLabel.Text = 'Fitting Style';
            app.FittingStyleDropDownLabel.FontColor = [11/255 71/255 111/255];

            % Create FittingStyleDropDown
            app.FittingStyleDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel,'FontColor',[11/255 71/255 111/255]);
            app.FittingStyleDropDown.BackgroundColor = [255/255 254/255 254/255];
            app.FittingStyleDropDown.Items = {'Concatenated', 'Separate'};
            app.FittingStyleDropDown.Enable = 'off';
            app.FittingStyleDropDown.Position = [161 145 129 22];
            app.FittingStyleDropDown.Value = 'Concatenated';

            % Create MRSFitRangeppmEditFieldLabel
            app.MRSFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.MRSFitRangeppmEditFieldLabel.HorizontalAlignment = 'right';
            app.MRSFitRangeppmEditFieldLabel.Position = [28 112 121 22];
            app.MRSFitRangeppmEditFieldLabel.Text = 'MRS Fit Range (ppm)';
            app.MRSFitRangeppmEditFieldLabel.FontColor = [11/255 71/255 111/255];

            % Create MRSFitRangeppmEditField
            app.MRSFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.MRSFitRangeppmEditField.Position = [161 112 129 22];
            app.MRSFitRangeppmEditField.Value = '0.2 4.2';

            % Create WaterFitRangeppmEditFieldLabel
            app.WaterFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.WaterFitRangeppmEditFieldLabel.HorizontalAlignment = 'right';
            app.WaterFitRangeppmEditFieldLabel.Position = [28 79 126 22];
            app.WaterFitRangeppmEditFieldLabel.Text = 'Water Fit Range (ppm)';
            app.WaterFitRangeppmEditFieldLabel.FontColor = [11/255 71/255 111/255];

            % Create WaterFitRangeppmEditField
            app.WaterFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.WaterFitRangeppmEditField.Position = [161 79 129 22];
            app.WaterFitRangeppmEditField.Value = '2.0 7.4';

            % Create BaselineknotspacingppmEditFieldLabel
            app.BaselineknotspacingppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.BaselineknotspacingppmEditFieldLabel.HorizontalAlignment = 'right';
            app.BaselineknotspacingppmEditFieldLabel.Position = [28 46 158 22];
            app.BaselineknotspacingppmEditFieldLabel.Text = 'Baseline knot spacing (ppm)';
            app.BaselineknotspacingppmEditFieldLabel.FontColor = [11/255 71/255 111/255];

            % Create BaselineknotspacingppmEditField
            app.BaselineknotspacingppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'numeric');
            app.BaselineknotspacingppmEditField.ValueDisplayFormat = '%11.1g';
            app.BaselineknotspacingppmEditField.HorizontalAlignment = 'left';
            app.BaselineknotspacingppmEditField.Position = [201 46 89 22];
            app.BaselineknotspacingppmEditField.Value = 0.4;

            % Create AddMMandLipbasisfunctionstofitCheckBox
            app.AddMMandLipbasisfunctionstofitCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.AddMMandLipbasisfunctionstofitCheckBox.Text = 'Add MM and Lip basis functions to fit';
            app.AddMMandLipbasisfunctionstofitCheckBox.Position = [28 14 225 22];
            app.AddMMandLipbasisfunctionstofitCheckBox.Value = true;
            app.AddMMandLipbasisfunctionstofitCheckBox.FontColor = [11/255 71/255 111/255];

            % Create SelectedMetabolitesPanel
            app.SelectedMetabolitesPanel = uipanel(app.UIFigure);
            app.SelectedMetabolitesPanel.ForegroundColor = [11/255 71/255 115/255];
            app.SelectedMetabolitesPanel.BackgroundColor = [255/255 254/255 254/255];
            app.SelectedMetabolitesPanel.TitlePosition = 'centertop';
            app.SelectedMetabolitesPanel.Title = 'Selected Metabolites';
            app.SelectedMetabolitesPanel.FontWeight = 'bold';
            app.SelectedMetabolitesPanel.FontSize = 16;
            app.SelectedMetabolitesPanel.Position = [5 10 707 150];

            % Create AlaCheckBox
            app.AlaCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AlaCheckBox.Enable = 'off';
            app.AlaCheckBox.Text = 'Ala';
            app.AlaCheckBox.Position = [15 90 39 22];
            app.AlaCheckBox.FontColor = [11/255 71/255 111/255];


            % Create AscCheckBox
            app.AscCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AscCheckBox.Enable = 'off';
            app.AscCheckBox.Text = 'Asc';
            app.AscCheckBox.Position = [75 90 42 22];
            app.AscCheckBox.Value = true;
            app.AscCheckBox.FontColor = [11/255 71/255 111/255];

            
            % Create AspCheckBox
            app.AspCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AspCheckBox.Enable = 'off';
            app.AspCheckBox.Text = 'Asp';
            app.AspCheckBox.Position = [135 90 43 22];
            app.AspCheckBox.Value = true;
            app.AspCheckBox.FontColor = [11/255 71/255 111/255];

            % Create bHBCheckBox
            app.bHBCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHBCheckBox.Enable = 'off';
            app.bHBCheckBox.Text = 'bHB';
            app.bHBCheckBox.Position = [195 90 46 22];
            app.bHBCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create bHGCheckBox
            app.bHGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHGCheckBox.Enable = 'off';
            app.bHGCheckBox.Text = 'bHG';
            app.bHGCheckBox.Position = [255 90 47 22];
            app.bHGCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create CitCheckBox
            app.CitCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CitCheckBox.Enable = 'off';
            app.CitCheckBox.Text = 'Cit';
            app.CitCheckBox.Position = [315 90 37 22];
            app.CitCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create CrCheckBox
            app.CrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCheckBox.Enable = 'off';
            app.CrCheckBox.Text = 'Cr';
            app.CrCheckBox.Position = [375 90 35 22];
            app.CrCheckBox.Value = true;
            app.CrCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create CrCH2CheckBox
            app.CrCH2CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCH2CheckBox.Enable = 'off';
            app.CrCH2CheckBox.Text = 'CrCH2';
            app.CrCH2CheckBox.Position = [425 90 59 22];
            app.CrCH2CheckBox.Value = true;
            app.CrCH2CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create EtOHCheckBox
            app.EtOHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.EtOHCheckBox.Enable = 'off';
            app.EtOHCheckBox.Text = 'EtOH';
            app.EtOHCheckBox.Position = [485 90 51 22];
            app.EtOHCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GABACheckBox
            app.GABACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GABACheckBox.Enable = 'off';
            app.GABACheckBox.Text = 'GABA';
            app.GABACheckBox.Position = [545 90 55 22];
            app.GABACheckBox.Value = true;
            app.GABACheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GPCCheckBox
            app.GPCCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GPCCheckBox.Enable = 'off';
            app.GPCCheckBox.Text = 'GPC';
            app.GPCCheckBox.Position = [605 90 48 22];
            app.GPCCheckBox.Value = true;
            app.GPCCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GSHCheckBox
            app.GSHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GSHCheckBox.Enable = 'off';
            app.GSHCheckBox.Text = 'GSH';
            app.GSHCheckBox.Position = [655 90 48 22];
            app.GSHCheckBox.Value = true;
            app.GSHCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GlcCheckBox
            app.GlcCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlcCheckBox.Enable = 'off';
            app.GlcCheckBox.Text = 'Glc';
            app.GlcCheckBox.Position = [15 65 40 22];
            app.GlcCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GlnCheckBox
            app.GlnCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlnCheckBox.Enable = 'off';
            app.GlnCheckBox.Text = 'Gln';
            app.GlnCheckBox.Position = [75 65 40 22];
            app.GlnCheckBox.Value = true;
            app.GlnCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GluCheckBox
            app.GluCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GluCheckBox.Enable = 'off';
            app.GluCheckBox.Text = 'Glu';
            app.GluCheckBox.Position = [135 65 40 22];
            app.GluCheckBox.Value = true;
            app.GluCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create GlyCheckBox
            app.GlyCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlyCheckBox.Enable = 'off';
            app.GlyCheckBox.Text = 'Gly';
            app.GlyCheckBox.Position = [195 65 40 22];
            app.GlyCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create H2OCheckBox
            app.H2OCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.H2OCheckBox.Enable = 'off';
            app.H2OCheckBox.Text = 'H2O';
            app.H2OCheckBox.Position = [255 65 46 22];
            app.H2OCheckBox.Value = true;
            app.H2OCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create InsCheckBox
            app.InsCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.InsCheckBox.Enable = 'off';
            app.InsCheckBox.Text = 'Ins';
            app.InsCheckBox.Position = [315 65 38 22];
            app.InsCheckBox.Value = true;
            app.InsCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create LacCheckBox
            app.LacCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.LacCheckBox.Enable = 'off';
            app.LacCheckBox.Text = 'Lac';
            app.LacCheckBox.Position = [375 65 42 22];
            app.LacCheckBox.Value = true;
            app.LacCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create NAACheckBox
            app.NAACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAACheckBox.Enable = 'off';
            app.NAACheckBox.Text = 'NAA';
            app.NAACheckBox.Position = [425 65 46 22];
            app.NAACheckBox.Value = true;
            app.NAACheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create NAAGCheckBox
            app.NAAGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAAGCheckBox.Enable = 'off';
            app.NAAGCheckBox.Text = 'NAAG';
            app.NAAGCheckBox.Position = [485 65 55 22];
            app.NAAGCheckBox.Value = true;
            app.NAAGCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create PChCheckBox
            app.PChCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PChCheckBox.Enable = 'off';
            app.PChCheckBox.Text = 'PCh';
            app.PChCheckBox.Position = [545 65 45 22];
            app.PChCheckBox.Value = true;
            app.PChCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create PCrCheckBox
            app.PCrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PCrCheckBox.Enable = 'off';
            app.PCrCheckBox.Text = 'PCr';
            app.PCrCheckBox.Position = [605 65 42 22];
            app.PCrCheckBox.Value = true;
            app.PCrCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create PECheckBox
            app.PECheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PECheckBox.Enable = 'off';
            app.PECheckBox.Text = 'PE';
            app.PECheckBox.Position = [665 65 37 22];
            app.PECheckBox.Value = true;
            app.PECheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create PhenylCheckBox
            app.PhenylCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PhenylCheckBox.Enable = 'off';
            app.PhenylCheckBox.Text = 'Phenyl';
            app.PhenylCheckBox.Position = [15 40 58 22];
            app.PhenylCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create ScylloCheckBox
            app.ScylloCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.ScylloCheckBox.Enable = 'off';
            app.ScylloCheckBox.Text = 'Scyllo';
            app.ScylloCheckBox.Position = [75 40 54 22];
            app.ScylloCheckBox.Value = true;
            app.ScylloCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create SerCheckBox
            app.SerCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.SerCheckBox.Enable = 'off';
            app.SerCheckBox.Text = 'Ser';
            app.SerCheckBox.Position = [135 40 40 22];
            app.SerCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create TauCheckBox
            app.TauCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TauCheckBox.Enable = 'off';
            app.TauCheckBox.Text = 'Tau';
            app.TauCheckBox.Position = [195 40 41 22];
            app.TauCheckBox.Value = true;
            app.TauCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create TyrosCheckBox
            app.TyrosCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TyrosCheckBox.Enable = 'off';
            app.TyrosCheckBox.Text = 'Tyros';
            app.TyrosCheckBox.Position = [255 40 50 22];
            app.TyrosCheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create MM09CheckBox
            app.MM09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM09CheckBox.Enable = 'off';
            app.MM09CheckBox.Text = 'MM09';
            app.MM09CheckBox.Position = [15 15 56 22];
            app.MM09CheckBox.Value = true;
            app.MM09CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create MM12CheckBox
            app.MM12CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM12CheckBox.Enable = 'off';
            app.MM12CheckBox.Text = 'MM12';
            app.MM12CheckBox.Position = [75 15 56 22];
            app.MM12CheckBox.Value = true;
            app.MM12CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create MM14CheckBox
            app.MM14CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM14CheckBox.Enable = 'off';
            app.MM14CheckBox.Text = 'MM14';
            app.MM14CheckBox.Position = [135 15 56 22];
            app.MM14CheckBox.Value = true;
            app.MM14CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create MM17CheckBox
            app.MM17CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM17CheckBox.Enable = 'off';
            app.MM17CheckBox.Text = 'MM17';
            app.MM17CheckBox.Position = [195 15 56 22];
            app.MM17CheckBox.Value = true;
            app.MM17CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create MM20CheckBox
            app.MM20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM20CheckBox.Enable = 'off';
            app.MM20CheckBox.Text = 'MM20';
            app.MM20CheckBox.Position = [255 15 56 22];
            app.MM20CheckBox.Value = true;
            app.MM20CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create Lip09CheckBox
            app.Lip09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip09CheckBox.Enable = 'off';
            app.Lip09CheckBox.Text = 'Lip09';
            app.Lip09CheckBox.Position = [315 15 52 22];
            app.Lip09CheckBox.Value = true;
            app.Lip09CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create Lip13CheckBox
            app.Lip13CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip13CheckBox.Enable = 'off';
            app.Lip13CheckBox.Text = 'Lip13';
            app.Lip13CheckBox.Position = [375 15 52 22];
            app.Lip13CheckBox.Value = true;
            app.Lip13CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create Lip20CheckBox
            app.Lip20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip20CheckBox.Enable = 'off';
            app.Lip20CheckBox.Text = 'Lip20';
            app.Lip20CheckBox.Position = [435 15 52 22];
            app.Lip20CheckBox.Value = true;
            app.Lip20CheckBox.FontColor = [11/255 71/255 111/255];
            
            % Create SpecifyMRSandStructuralImagingFilesPanel
            app.SpecifyMRSandStructuralImagingFilesPanel = uipanel(app.UIFigure);
            app.SpecifyMRSandStructuralImagingFilesPanel.ForegroundColor = [11/255 71/255 115/255];
            app.SpecifyMRSandStructuralImagingFilesPanel.BackgroundColor = [255/255 254/255 254/255];
            app.SpecifyMRSandStructuralImagingFilesPanel.Title = '3. Specify MRS and Structural Imaging Files';
            app.SpecifyMRSandStructuralImagingFilesPanel.FontWeight = 'bold';
            app.SpecifyMRSandStructuralImagingFilesPanel.FontSize = 16;
            app.SpecifyMRSandStructuralImagingFilesPanel.Position = [362 345 350 249];

            % Create MRSDataButton
            app.MRSDataButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.MRSDataButton.ButtonPushedFcn = createCallbackFcn(app, @MRSDataButtonPushed, true);
            app.MRSDataButton.BackgroundColor = [255/255 254/255 254/255];
            app.MRSDataButton.FontWeight = 'bold';
            app.MRSDataButton.Position = [10 151 118 23];
            app.MRSDataButton.Text = 'metabolite spectra';
            app.MRSDataButton.FontColor = [11/255 71/255 111/255];

            % Create H2OReferenceButton
            app.H2OReferenceButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.H2OReferenceButton.ButtonPushedFcn = createCallbackFcn(app, @H2OReferenceButtonPushed, true);
            app.H2OReferenceButton.BackgroundColor = [255/255 254/255 254/255];
            app.H2OReferenceButton.FontWeight = 'bold';
            app.H2OReferenceButton.FontColor = [11/255 71/255 111/255];
            app.H2OReferenceButton.Enable = 'off';
            app.H2OReferenceButton.Position = [10 119 118 23];
            app.H2OReferenceButton.Text = 'water reference';

            % Create H2OShortTEButton
            app.H2OShortTEButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.H2OShortTEButton.ButtonPushedFcn = createCallbackFcn(app, @H2OShortTEButtonPushed, true);
            app.H2OShortTEButton.BackgroundColor = [255/255 254/255 254/255];
            app.H2OShortTEButton.FontWeight = 'bold';
            app.H2OShortTEButton.FontColor = [11/255 71/255 111/255];
            app.H2OShortTEButton.Enable = 'off';
            app.H2OShortTEButton.Position = [10 84 118 23];
            app.H2OShortTEButton.Text = 'short TE water';

            % Create MetaboliteNulledButton
            app.MetaboliteNulledButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.MetaboliteNulledButton.ButtonPushedFcn = createCallbackFcn(app, @MetaboliteNulledButtonPushed, true);
            app.MetaboliteNulledButton.BackgroundColor = [255/255 254/255 254/255];
            app.MetaboliteNulledButton.FontWeight = 'bold';
            app.MetaboliteNulledButton.FontColor = [11/255 71/255 111/255];
            app.MetaboliteNulledButton.Enable = 'off';
            app.MetaboliteNulledButton.Position = [10 52 118 23];
            app.MetaboliteNulledButton.Text = 'metabolite-nulled';

            % Create T1DataniftiniiButton
            app.T1DataniftiniiButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.T1DataniftiniiButton.ButtonPushedFcn = createCallbackFcn(app, @T1DataniftiniiButtonPushed, true);
            app.T1DataniftiniiButton.BackgroundColor = [255/255 254/255 254/255];
            app.T1DataniftiniiButton.FontWeight = 'bold';
            app.T1DataniftiniiButton.FontColor = [11/255 71/255 111/255];
            app.T1DataniftiniiButton.Enable = 'off';
            app.T1DataniftiniiButton.Position = [11 18 117 23];
            app.T1DataniftiniiButton.Text = 'T1 anatomical';

            % Create NumberofdatasetsEditFieldLabel
            app.NumberofdatasetsEditFieldLabel = uilabel(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.NumberofdatasetsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofdatasetsEditFieldLabel.Position = [18 191 112 22];
            app.NumberofdatasetsEditFieldLabel.Text = 'Number of datasets';
            app.NumberofdatasetsEditFieldLabel.FontColor = [11/255 71/255 111/255];

            % Create NumberofdatasetsEditField
            app.NumberofdatasetsEditField = uieditfield(app.SpecifyMRSandStructuralImagingFilesPanel, 'numeric');
            app.NumberofdatasetsEditField.Limits = [1 Inf];
            app.NumberofdatasetsEditField.ValueDisplayFormat = '%.0f';
            app.NumberofdatasetsEditField.HorizontalAlignment = 'center';
            app.NumberofdatasetsEditField.Position = [145 191 100 22];
            app.NumberofdatasetsEditField.Value = 1;

            % Create MRSDataText
            app.MRSDataText = uitextarea(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.MRSDataText.Editable = 'off';
            app.MRSDataText.Position = [144 150 188 25];

            % Create H2OShortTEText
            app.H2OShortTEText = uitextarea(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.H2OShortTEText.Editable = 'off';
            app.H2OShortTEText.Position = [144 83 188 25];

            % Create H2OReferenceText
            app.H2OReferenceText = uitextarea(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.H2OReferenceText.Editable = 'off';
            app.H2OReferenceText.Position = [144 118 188 25];

            % Create MetaboliteNulledText
            app.MetaboliteNulledText = uitextarea(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.MetaboliteNulledText.Editable = 'off';
            app.MetaboliteNulledText.Position = [144 51 188 25];

            % Create T1DataText
            app.T1DataText = uitextarea(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.T1DataText.Editable = 'off';
            app.T1DataText.Position = [144 17 188 25];

            % Create SpecifyOutputFolderPanel
            app.SpecifyOutputFolderPanel = uipanel(app.UIFigure);
            app.SpecifyOutputFolderPanel.ForegroundColor = [11/255 71/255 115/255];
            app.SpecifyOutputFolderPanel.BackgroundColor = [255/255 254/255 254/255];
            app.SpecifyOutputFolderPanel.Title = '4. Specify Output Folder';
            app.SpecifyOutputFolderPanel.FontWeight = 'bold';
            app.SpecifyOutputFolderPanel.FontSize = 16;
            app.SpecifyOutputFolderPanel.Position = [362 225 349 107];

            % Create OutputFolderEditField
            app.OutputFolderEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.OutputFolderEditField.Editable = 'off';
            app.OutputFolderEditField.Position = [144 48 183 22];

            % Create OutputFolderButton
            app.OutputFolderButton = uibutton(app.SpecifyOutputFolderPanel, 'push');
            app.OutputFolderButton.ButtonPushedFcn = createCallbackFcn(app, @OutputFolderButtonPushed, true);
            app.OutputFolderButton.BackgroundColor = [255/255 254/255 254/255];
            app.OutputFolderButton.FontWeight = 'bold';
            app.OutputFolderButton.FontColor = [11/255 71/255 111/255];
            app.OutputFolderButton.Position = [9 47 119 23];
            app.OutputFolderButton.Text = 'Output Folder';

            % Create JobNameEditFieldLabel
            app.JobNameEditFieldLabel = uilabel(app.SpecifyOutputFolderPanel);
            app.JobNameEditFieldLabel.HorizontalAlignment = 'right';
            app.JobNameEditFieldLabel.Position = [18 12 61 22];
            app.JobNameEditFieldLabel.Text = 'Job Name';
            app.JobNameEditFieldLabel.FontColor = [11/255 71/255 111/255];

            % Create JobNameEditField
            app.JobNameEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.JobNameEditField.Position = [144 12 183 22];

            % Create CANCELButton
            app.CANCELButton = uibutton(app.UIFigure, 'push');
            app.CANCELButton.ButtonPushedFcn = createCallbackFcn(app, @CANCELButtonPushed, true);
            app.CANCELButton.BackgroundColor = [255/255 254/255 254/255];
            app.CANCELButton.FontSize = 16;
            app.CANCELButton.FontWeight = 'bold';
            app.CANCELButton.FontColor = [11/255 71/255 111/255];
            app.CANCELButton.Position = [550 170 160 43];           
            app.CANCELButton.Text = 'Cancel';

            % Create CREATEJOBButton
            app.CREATEJOBButton = uibutton(app.UIFigure, 'push');
            app.CREATEJOBButton.ButtonPushedFcn = createCallbackFcn(app, @CREATEJOBButtonPushed, true);
            app.CREATEJOBButton.BackgroundColor = [255/255 254/255 254/255];
            app.CREATEJOBButton.FontSize = 16;
            app.CREATEJOBButton.FontWeight = 'bold';
            app.CREATEJOBButton.FontColor = [11/255 71/255 111/255];
            app.CREATEJOBButton.Position = [367 170 160 43];
            app.CREATEJOBButton.Text = 'Create JobFile';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CreateOspreyJob_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end