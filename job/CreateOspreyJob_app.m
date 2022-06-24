classdef CreateOspreyJob_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        InteractiveOspreyjobfilegeneratorUIFigure  matlab.ui.Figure
        FWHMMM3coEditField              matlab.ui.control.NumericEditField
        FWHMMM3coLabel                  matlab.ui.control.Label
        basissetfileButton              matlab.ui.control.Button
        BasisSetEditField               matlab.ui.control.EditField
        StatcsvFileButton               matlab.ui.control.Button
        Image                           matlab.ui.control.Image
        CREATEJOBButton                 matlab.ui.control.Button
        CANCELButton                    matlab.ui.control.Button
        StatcsvEditField                matlab.ui.control.EditField
        SpecifyOutputFolderPanel        matlab.ui.container.Panel
        JobNameEditField                matlab.ui.control.EditField
        JobNameEditFieldLabel           matlab.ui.control.Label
        OutputFolderButton              matlab.ui.control.Button
        OutputFolderEditField           matlab.ui.control.EditField
        SpecifyMRSandAnatomicalImagingFilesPanel  matlab.ui.container.Panel
        T1DataText                      matlab.ui.control.TextArea
        MetaboliteNulledText            matlab.ui.control.TextArea
        H2OReferenceText                matlab.ui.control.TextArea
        H2OShortTEText                  matlab.ui.control.TextArea
        MRSDataText                     matlab.ui.control.TextArea
        NumberofdatasetsEditField       matlab.ui.control.NumericEditField
        NumberofdatasetsEditFieldLabel  matlab.ui.control.Label
        T1DataniftiniiButton            matlab.ui.control.Button
        MetaboliteNulledButton          matlab.ui.control.Button
        H2OShortTEButton                matlab.ui.control.Button
        H2OReferenceButton              matlab.ui.control.Button
        MRSDataButton                   matlab.ui.control.Button
        SelectedMetabolitesPanel        matlab.ui.container.Panel
        CystatCheckBox                  matlab.ui.control.CheckBox
        MM_CSOCheckBox                  matlab.ui.control.CheckBox
        MM_PCCCheckBox                  matlab.ui.control.CheckBox
        MMexpCheckBox                   matlab.ui.control.CheckBox
        MM42CheckBox                    matlab.ui.control.CheckBox
        MM40CheckBox                    matlab.ui.control.CheckBox
        MM38CheckBox                    matlab.ui.control.CheckBox
        MM37CheckBox                    matlab.ui.control.CheckBox
        Lip20CheckBox                   matlab.ui.control.CheckBox
        Lip13CheckBox                   matlab.ui.control.CheckBox
        Lip09CheckBox                   matlab.ui.control.CheckBox
        MM20CheckBox                    matlab.ui.control.CheckBox
        MM17CheckBox                    matlab.ui.control.CheckBox
        MM14CheckBox                    matlab.ui.control.CheckBox
        MM12CheckBox                    matlab.ui.control.CheckBox
        MM09CheckBox                    matlab.ui.control.CheckBox
        TyrosCheckBox                   matlab.ui.control.CheckBox
        NAAGCheckBox                    matlab.ui.control.CheckBox
        TauCheckBox                     matlab.ui.control.CheckBox
        SerCheckBox                     matlab.ui.control.CheckBox
        ScylloCheckBox                  matlab.ui.control.CheckBox
        PhenylCheckBox                  matlab.ui.control.CheckBox
        PECheckBox                      matlab.ui.control.CheckBox
        PCrCheckBox                     matlab.ui.control.CheckBox
        PChCheckBox                     matlab.ui.control.CheckBox
        GluCheckBox                     matlab.ui.control.CheckBox
        NAACheckBox                     matlab.ui.control.CheckBox
        LacCheckBox                     matlab.ui.control.CheckBox
        InsCheckBox                     matlab.ui.control.CheckBox
        H2OCheckBox                     matlab.ui.control.CheckBox
        GlyCheckBox                     matlab.ui.control.CheckBox
        GlnCheckBox                     matlab.ui.control.CheckBox
        GlcCheckBox                     matlab.ui.control.CheckBox
        GSHCheckBox                     matlab.ui.control.CheckBox
        GPCCheckBox                     matlab.ui.control.CheckBox
        GABACheckBox                    matlab.ui.control.CheckBox
        EtOHCheckBox                    matlab.ui.control.CheckBox
        CrCH2CheckBox                   matlab.ui.control.CheckBox
        CrCheckBox                      matlab.ui.control.CheckBox
        CitCheckBox                     matlab.ui.control.CheckBox
        bHGCheckBox                     matlab.ui.control.CheckBox
        bHBCheckBox                     matlab.ui.control.CheckBox
        AspCheckBox                     matlab.ui.control.CheckBox
        AscCheckBox                     matlab.ui.control.CheckBox
        AlaCheckBox                     matlab.ui.control.CheckBox
        SpecifyDataHandlingandModelingOptionsPanel  matlab.ui.container.Panel
        SavePDFCheckBox                 matlab.ui.control.CheckBox
        unstablewaterCheckBox           matlab.ui.control.CheckBox
        ECCmmMRSCheckBox                matlab.ui.control.CheckBox
        ECCmetaboliteMRSCheckBox        matlab.ui.control.CheckBox
        SaveNIfTiCheckBox               matlab.ui.control.CheckBox
        SubspectraAligmentDropDown      matlab.ui.control.DropDown
        SubspectraAligmentDropDownLabel  matlab.ui.control.Label
        SpectralRegistrationDropDown    matlab.ui.control.DropDown
        SpectralRegistrationDropDownLabel  matlab.ui.control.Label
        AddMMandLipbasisfunctionstofitCheckBox  matlab.ui.control.CheckBox
        BaselineknotspacingppmEditField  matlab.ui.control.NumericEditField
        BaselineknotspacingppmEditFieldLabel  matlab.ui.control.Label
        WaterFitRangeppmEditField       matlab.ui.control.EditField
        WaterFitRangeppmEditFieldLabel  matlab.ui.control.Label
        MRSFitRangeppmEditField         matlab.ui.control.EditField
        MRSFitRangeppmEditFieldLabel    matlab.ui.control.Label
        FittingStyleDropDown            matlab.ui.control.DropDown
        FittingStyleDropDownLabel       matlab.ui.control.Label
        IncludedMetabolitesDropDown     matlab.ui.control.DropDown
        IncludedMetabolitesDropDownLabel  matlab.ui.control.Label
        FittingAlgorithmDropDown        matlab.ui.control.DropDown
        FittingAlgorithmDropDownLabel   matlab.ui.control.Label
        SaveVendorCheckBox              matlab.ui.control.CheckBox
        SaveJMRUICheckBox               matlab.ui.control.CheckBox
        SaveLCMCheckBox                 matlab.ui.control.CheckBox
        SpecifySequenceInformationPanel  matlab.ui.container.Panel
        MM3coDropDown                   matlab.ui.control.DropDown
        MM3coDropDownLabel              matlab.ui.control.Label
        DataScenarioDropDown            matlab.ui.control.DropDown
        DataScenarioDropDownLabel       matlab.ui.control.Label
        EditingTargetsDropDown          matlab.ui.control.DropDown
        EditingTargetsDropDownLabel     matlab.ui.control.Label
        SequenceTypeDropDown            matlab.ui.control.DropDown
        SequenceTypeDropDownLabel       matlab.ui.control.Label
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
                    app.FWHMMM3coEditField.Enable = 'Off';
                    app.MM3coDropDown.Enable = 'Off';
                    app.MM3coDropDownLabel.Enable = 'Off';
                case 'MEGA'
                    app.EditingTargetsDropDown.Items = {'GABA','GSH','Lac','PE322','PE398'};
                    app.EditingTargetsDropDown.Value = 'GABA';
                    
                    app.FittingStyleDropDown.Enable = 'Off';
                    app.FittingStyleDropDownLabel.Enable = 'Off';
                    app.FWHMMM3coEditField.Enable = 'On';
                    app.MM3coDropDown.Enable = 'On';
                    app.MM3coDropDownLabel.Enable = 'On';
                case 'HERMES'
                    app.EditingTargetsDropDown.Items = {'GABA, GSH'};
                    app.EditingTargetsDropDown.Value = 'GABA, GSH';
                    
                    app.FittingStyleDropDown.Enable = 'Off';
                    app.FittingStyleDropDownLabel.Enable = 'Off';
                    app.FWHMMM3coEditField.Enable = 'On';
                    app.MM3coDropDown.Enable = 'On';
                    app.MM3coDropDownLabel.Enable = 'On';
                case 'HERCULES'
                    app.EditingTargetsDropDown.Items = {'GABA, GSH'};
                    app.EditingTargetsDropDown.Value = 'GABA, GSH';
                    
                    app.FittingStyleDropDown.Enable = 'Off';
                    app.FittingStyleDropDownLabel.Enable = 'Off';
                    app.FWHMMM3coEditField.Enable = 'On';
                    app.MM3coDropDown.Enable = 'On';
                    app.MM3coDropDownLabel.Enable = 'On';
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
            info = 'Please select the MRS file to read';
            
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
        end

        % Button pushed function: H2OReferenceButton
        function H2OReferenceButtonPushed(app, event)
            info = 'Please select the water reference file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            h2oreffiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2oreffiles(i,:)};
            end
            
            app.H2OReferenceText.Value = filelist;
        end

        % Button pushed function: H2OShortTEButton
        function H2OShortTEButtonPushed(app, event)
            info = 'Please select the water short-TE file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            h2ostefiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2ostefiles(i,:)};
            end
            
            app.H2OShortTEText.Value = filelist;
        end

        % Button pushed function: MetaboliteNulledButton
        function MetaboliteNulledButtonPushed(app, event)
            info = 'Please select the metabolite-nulled file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            metnulfiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} metnulfiles(i,:)};
            end
            
            app.MetaboliteNulledText.Value = filelist;
        end

        % Button pushed function: T1DataniftiniiButton
        function T1DataniftiniiButtonPushed(app, event)
            info = 'Please select the T1 Data file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            t1imfiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} t1imfiles(i,:)};
            end
            
            app.T1DataText.Value = filelist;
        end

        % Button pushed function: OutputFolderButton
        function OutputFolderButtonPushed(app, event)
            info = 'Please select the output folder';
            pathname=uigetdir('*.*',info);
            
            app.OutputFolderEditField.Value = pathname;
        end

        % Button pushed function: CANCELButton
        function CANCELButtonPushed(app, event)
            delete(app.InteractiveOspreyjobfilegeneratorUIFigure)
        end

        % Button pushed function: CREATEJOBButton
        function CREATEJOBButtonPushed(app, event)
            jobm = osp_create_job_file(app);
            
            delete(app.InteractiveOspreyjobfilegeneratorUIFigure)
            
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

        % Button pushed function: StatcsvFileButton
        function StatcsvFileButtonPushed(app, event)
            info = 'Please select stat csv file to read';
            
            ndata = 1;
            
            csvfiles = spm_select(ndata,'any',info,{},pwd,'.csv','1');
            
            app.StatcsvEditField.Value = csvfiles(1,:);
        end

        % Button pushed function: basissetfileButton
        function basissetfileButtonPushed(app, event)
            info = 'Select a .mat basis file to overwrite the automatic basis set selection';
            
            ndata = 1;
            
            basisfiles = spm_select(ndata,'any',info,{},pwd,'.mat','1');
            
            app.BasisSetEditField.Value = basisfiles(1,:);
        
        end

        % Value changed function: EditingTargetsDropDown
        function EditingTargetsDropDownValueChanged(app, event)
            value = app.EditingTargetsDropDown.Value;
            switch value
                case 'GABA'
                    app.FWHMMM3coEditField.Enable = 'On';
                    app.MM3coDropDown.Enable = 'On';
                    app.MM3coDropDownLabel.Enable = 'On';
                case 'GABA, GSH'
                    app.FWHMMM3coEditField.Enable = 'On';
                    app.MM3coDropDown.Enable = 'On';
                    app.MM3coDropDownLabel.Enable = 'On';
                otherwise
                    app.FWHMMM3coEditField.Enable = 'Off';
                    app.MM3coDropDown.Enable = 'Off';
                    app.MM3coDropDownLabel.Enable = 'Off';

            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create InteractiveOspreyjobfilegeneratorUIFigure and hide until all components are created
            app.InteractiveOspreyjobfilegeneratorUIFigure = uifigure('Visible', 'off');
            app.InteractiveOspreyjobfilegeneratorUIFigure.Color = [1 0.9882 0.9882];
            app.InteractiveOspreyjobfilegeneratorUIFigure.Position = [100 100 643 850];
            app.InteractiveOspreyjobfilegeneratorUIFigure.Name = 'Interactive Osprey jobfile generator';

            % Create SpecifySequenceInformationPanel
            app.SpecifySequenceInformationPanel = uipanel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.SpecifySequenceInformationPanel.Tooltip = {'Specify Sequence Information here.'};
            app.SpecifySequenceInformationPanel.ForegroundColor = [0.0392 0.2706 0.4314];
            app.SpecifySequenceInformationPanel.BorderType = 'none';
            app.SpecifySequenceInformationPanel.Title = '1. Specify Sequence Information';
            app.SpecifySequenceInformationPanel.BackgroundColor = [1 0.9882 0.9882];
            app.SpecifySequenceInformationPanel.FontWeight = 'bold';
            app.SpecifySequenceInformationPanel.FontSize = 15;
            app.SpecifySequenceInformationPanel.Position = [6 748 630 95];

            % Create SequenceTypeDropDownLabel
            app.SequenceTypeDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.SequenceTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SequenceTypeDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.SequenceTypeDropDownLabel.Position = [9 39 88 22];
            app.SequenceTypeDropDownLabel.Text = 'Sequence Type';

            % Create SequenceTypeDropDown
            app.SequenceTypeDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.SequenceTypeDropDown.Items = {'unedited', 'MEGA', 'HERMES', 'HERCULES'};
            app.SequenceTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SequenceTypeDropDownValueChanged, true);
            app.SequenceTypeDropDown.Tooltip = {'Select the sequence type. The basis set will be determined automatically.'};
            app.SequenceTypeDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.SequenceTypeDropDown.Position = [112 39 129 22];
            app.SequenceTypeDropDown.Value = 'unedited';

            % Create EditingTargetsDropDownLabel
            app.EditingTargetsDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.EditingTargetsDropDownLabel.HorizontalAlignment = 'right';
            app.EditingTargetsDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.EditingTargetsDropDownLabel.Position = [252 39 85 22];
            app.EditingTargetsDropDownLabel.Text = 'Editing Targets';

            % Create EditingTargetsDropDown
            app.EditingTargetsDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.EditingTargetsDropDown.Items = {'none'};
            app.EditingTargetsDropDown.ValueChangedFcn = createCallbackFcn(app, @EditingTargetsDropDownValueChanged, true);
            app.EditingTargetsDropDown.Tooltip = {'Select the editing target. If yours isn''t availabe please contact the developers.'};
            app.EditingTargetsDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.EditingTargetsDropDown.Position = [355 39 129 22];
            app.EditingTargetsDropDown.Value = 'none';

            % Create DataScenarioDropDownLabel
            app.DataScenarioDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.DataScenarioDropDownLabel.HorizontalAlignment = 'right';
            app.DataScenarioDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.DataScenarioDropDownLabel.Position = [16 7 81 22];
            app.DataScenarioDropDownLabel.Text = 'Data Scenario';

            % Create DataScenarioDropDown
            app.DataScenarioDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.DataScenarioDropDown.Items = {'invivo', 'phantom', 'PRIAM'};
            app.DataScenarioDropDown.Tooltip = {'Select a data scenario.'};
            app.DataScenarioDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.DataScenarioDropDown.Position = [112 7 129 22];
            app.DataScenarioDropDown.Value = 'invivo';

            % Create MM3coDropDownLabel
            app.MM3coDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.MM3coDropDownLabel.HorizontalAlignment = 'right';
            app.MM3coDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.MM3coDropDownLabel.Position = [254 7 83 22];
            app.MM3coDropDownLabel.Text = 'MM3co model';

            % Create MM3coDropDown
            app.MM3coDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.MM3coDropDown.Items = {'3to2MM', '3to2MMsoft', '1to1GABA', '1to1GABAsoft', 'freeGauss', 'fixedGauss', 'none'};
            app.MM3coDropDown.Enable = 'off';
            app.MM3coDropDown.Tooltip = {'Select a model for the co-edited MMs in GABA-edited difference spetra.'};
            app.MM3coDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.MM3coDropDown.Position = [355 7 129 22];
            app.MM3coDropDown.Value = '3to2MM';

            % Create SpecifyDataHandlingandModelingOptionsPanel
            app.SpecifyDataHandlingandModelingOptionsPanel = uipanel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.SpecifyDataHandlingandModelingOptionsPanel.Tooltip = {'Specify Preprocessing and Modeling options here.'};
            app.SpecifyDataHandlingandModelingOptionsPanel.ForegroundColor = [0.0392 0.2706 0.4314];
            app.SpecifyDataHandlingandModelingOptionsPanel.BorderType = 'none';
            app.SpecifyDataHandlingandModelingOptionsPanel.Title = '2. Specify Data Handling and Modeling Options';
            app.SpecifyDataHandlingandModelingOptionsPanel.BackgroundColor = [1 0.9882 0.9882];
            app.SpecifyDataHandlingandModelingOptionsPanel.FontWeight = 'bold';
            app.SpecifyDataHandlingandModelingOptionsPanel.FontSize = 15;
            app.SpecifyDataHandlingandModelingOptionsPanel.Position = [6 482 627 261];

            % Create SaveLCMCheckBox
            app.SaveLCMCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveLCMCheckBox.Tooltip = {'Save LCM .RAW files?'};
            app.SaveLCMCheckBox.Text = 'Save LCM';
            app.SaveLCMCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SaveLCMCheckBox.Position = [459 95 78 22];

            % Create SaveJMRUICheckBox
            app.SaveJMRUICheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveJMRUICheckBox.Tooltip = {'Save jMRUI .txt files?'};
            app.SaveJMRUICheckBox.Text = 'Save JMRUI';
            app.SaveJMRUICheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SaveJMRUICheckBox.Position = [459 43 89 22];

            % Create SaveVendorCheckBox
            app.SaveVendorCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveVendorCheckBox.Tooltip = {'Save vendor native files?'};
            app.SaveVendorCheckBox.Text = 'Save Vendor';
            app.SaveVendorCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SaveVendorCheckBox.Position = [459 69 90 22];

            % Create FittingAlgorithmDropDownLabel
            app.FittingAlgorithmDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingAlgorithmDropDownLabel.HorizontalAlignment = 'right';
            app.FittingAlgorithmDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.FittingAlgorithmDropDownLabel.Position = [9 208 94 22];
            app.FittingAlgorithmDropDownLabel.Text = 'Fitting Algorithm';

            % Create FittingAlgorithmDropDown
            app.FittingAlgorithmDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingAlgorithmDropDown.Items = {'Osprey', 'LCModel'};
            app.FittingAlgorithmDropDown.Tooltip = {'Select the model algorithm.'};
            app.FittingAlgorithmDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.FittingAlgorithmDropDown.Position = [142 208 129 22];
            app.FittingAlgorithmDropDown.Value = 'Osprey';

            % Create IncludedMetabolitesDropDownLabel
            app.IncludedMetabolitesDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.IncludedMetabolitesDropDownLabel.HorizontalAlignment = 'right';
            app.IncludedMetabolitesDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.IncludedMetabolitesDropDownLabel.Position = [9 174 118 22];
            app.IncludedMetabolitesDropDownLabel.Text = 'Included Metabolites';

            % Create IncludedMetabolitesDropDown
            app.IncludedMetabolitesDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.IncludedMetabolitesDropDown.Items = {'default', 'full', 'custom'};
            app.IncludedMetabolitesDropDown.ValueChangedFcn = createCallbackFcn(app, @IncludedMetabolitesDropDownValueChanged, true);
            app.IncludedMetabolitesDropDown.Tooltip = {'Change the metabolite list for the LCM.'};
            app.IncludedMetabolitesDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.IncludedMetabolitesDropDown.Position = [142 174 129 22];
            app.IncludedMetabolitesDropDown.Value = 'default';

            % Create FittingStyleDropDownLabel
            app.FittingStyleDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingStyleDropDownLabel.HorizontalAlignment = 'right';
            app.FittingStyleDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.FittingStyleDropDownLabel.Enable = 'off';
            app.FittingStyleDropDownLabel.Position = [9 140 69 22];
            app.FittingStyleDropDownLabel.Text = 'Fitting Style';

            % Create FittingStyleDropDown
            app.FittingStyleDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingStyleDropDown.Items = {'Concatenated', 'Separate'};
            app.FittingStyleDropDown.Enable = 'off';
            app.FittingStyleDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.FittingStyleDropDown.Position = [142 140 129 22];
            app.FittingStyleDropDown.Value = 'Separate';

            % Create MRSFitRangeppmEditFieldLabel
            app.MRSFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.MRSFitRangeppmEditFieldLabel.HorizontalAlignment = 'right';
            app.MRSFitRangeppmEditFieldLabel.FontColor = [0.0392 0.2706 0.4314];
            app.MRSFitRangeppmEditFieldLabel.Position = [9 104 121 22];
            app.MRSFitRangeppmEditFieldLabel.Text = 'MRS Fit Range (ppm)';

            % Create MRSFitRangeppmEditField
            app.MRSFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.MRSFitRangeppmEditField.FontColor = [0.0392 0.2706 0.4314];
            app.MRSFitRangeppmEditField.Tooltip = {'Change the fit range for the metabolite data.'};
            app.MRSFitRangeppmEditField.Position = [142 104 129 22];
            app.MRSFitRangeppmEditField.Value = '0.5 4.0';

            % Create WaterFitRangeppmEditFieldLabel
            app.WaterFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.WaterFitRangeppmEditFieldLabel.HorizontalAlignment = 'right';
            app.WaterFitRangeppmEditFieldLabel.FontColor = [0.0392 0.2706 0.4314];
            app.WaterFitRangeppmEditFieldLabel.Position = [9 70 126 22];
            app.WaterFitRangeppmEditFieldLabel.Text = 'Water Fit Range (ppm)';

            % Create WaterFitRangeppmEditField
            app.WaterFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.WaterFitRangeppmEditField.FontColor = [0.0392 0.2706 0.4314];
            app.WaterFitRangeppmEditField.Tooltip = {'Change the fit range for the water reference scan.'};
            app.WaterFitRangeppmEditField.Position = [142 70 129 22];
            app.WaterFitRangeppmEditField.Value = '2.0 7.4';

            % Create BaselineknotspacingppmEditFieldLabel
            app.BaselineknotspacingppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.BaselineknotspacingppmEditFieldLabel.HorizontalAlignment = 'right';
            app.BaselineknotspacingppmEditFieldLabel.FontColor = [0.0392 0.2706 0.4314];
            app.BaselineknotspacingppmEditFieldLabel.Position = [9 39 158 22];
            app.BaselineknotspacingppmEditFieldLabel.Text = 'Baseline knot spacing (ppm)';

            % Create BaselineknotspacingppmEditField
            app.BaselineknotspacingppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'numeric');
            app.BaselineknotspacingppmEditField.Limits = [0 Inf];
            app.BaselineknotspacingppmEditField.ValueDisplayFormat = '%11.1g';
            app.BaselineknotspacingppmEditField.HorizontalAlignment = 'left';
            app.BaselineknotspacingppmEditField.FontColor = [0.0392 0.2706 0.4314];
            app.BaselineknotspacingppmEditField.Tooltip = {'Change the baseline knot spcaing of the spline functions.'};
            app.BaselineknotspacingppmEditField.Position = [182 39 89 22];
            app.BaselineknotspacingppmEditField.Value = 0.4;

            % Create AddMMandLipbasisfunctionstofitCheckBox
            app.AddMMandLipbasisfunctionstofitCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.AddMMandLipbasisfunctionstofitCheckBox.Text = 'Add MM and Lip basis functions to fit';
            app.AddMMandLipbasisfunctionstofitCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.AddMMandLipbasisfunctionstofitCheckBox.Position = [16 1 225 22];
            app.AddMMandLipbasisfunctionstofitCheckBox.Value = true;

            % Create SpectralRegistrationDropDownLabel
            app.SpectralRegistrationDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SpectralRegistrationDropDownLabel.HorizontalAlignment = 'right';
            app.SpectralRegistrationDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.SpectralRegistrationDropDownLabel.Position = [308 208 118 22];
            app.SpectralRegistrationDropDownLabel.Text = 'Spectral Registration';

            % Create SpectralRegistrationDropDown
            app.SpectralRegistrationDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SpectralRegistrationDropDown.Items = {'RobSpecReg', 'ProbSpecReg', 'RestrSpecReg', 'none'};
            app.SpectralRegistrationDropDown.Tooltip = {'Change the spectral registration algorithm for the alignemnt of the transients.'};
            app.SpectralRegistrationDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.SpectralRegistrationDropDown.Position = [456 208 129 22];
            app.SpectralRegistrationDropDown.Value = 'RobSpecReg';

            % Create SubspectraAligmentDropDownLabel
            app.SubspectraAligmentDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SubspectraAligmentDropDownLabel.HorizontalAlignment = 'right';
            app.SubspectraAligmentDropDownLabel.FontColor = [0.0392 0.2706 0.4314];
            app.SubspectraAligmentDropDownLabel.Position = [308 174 118 22];
            app.SubspectraAligmentDropDownLabel.Text = 'Subspectra Aligment';

            % Create SubspectraAligmentDropDown
            app.SubspectraAligmentDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SubspectraAligmentDropDown.Items = {'L2Norm', 'L1Norm', 'none'};
            app.SubspectraAligmentDropDown.Tooltip = {'Change the subspectra alignment algorithm.'};
            app.SubspectraAligmentDropDown.FontColor = [0.0392 0.2706 0.4314];
            app.SubspectraAligmentDropDown.Position = [456 174 129 22];
            app.SubspectraAligmentDropDown.Value = 'L2Norm';

            % Create SaveNIfTiCheckBox
            app.SaveNIfTiCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveNIfTiCheckBox.Tooltip = {'Save NIfTi MRS files?'};
            app.SaveNIfTiCheckBox.Text = 'Save NIfTi';
            app.SaveNIfTiCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SaveNIfTiCheckBox.Position = [552 69 77 22];

            % Create ECCmetaboliteMRSCheckBox
            app.ECCmetaboliteMRSCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.ECCmetaboliteMRSCheckBox.Tooltip = {'Do eddy current correction on metabolite data?'};
            app.ECCmetaboliteMRSCheckBox.Text = 'ECC metabolite MRS';
            app.ECCmetaboliteMRSCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.ECCmetaboliteMRSCheckBox.Position = [313 95 136 22];
            app.ECCmetaboliteMRSCheckBox.Value = true;

            % Create ECCmmMRSCheckBox
            app.ECCmmMRSCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.ECCmmMRSCheckBox.Tooltip = {'Do eddy current correction on metabolite-nulled data?'};
            app.ECCmmMRSCheckBox.Text = 'ECC mm MRS';
            app.ECCmmMRSCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.ECCmmMRSCheckBox.Position = [313 64 100 22];
            app.ECCmmMRSCheckBox.Value = true;

            % Create unstablewaterCheckBox
            app.unstablewaterCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.unstablewaterCheckBox.Tooltip = {'Do eddy current correction on metabolite data?'};
            app.unstablewaterCheckBox.Text = 'unstable water';
            app.unstablewaterCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.unstablewaterCheckBox.Position = [458 140 101 22];

            % Create SavePDFCheckBox
            app.SavePDFCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SavePDFCheckBox.Tooltip = {'Save LCM .RAW files?'};
            app.SavePDFCheckBox.Text = 'Save PDF';
            app.SavePDFCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SavePDFCheckBox.Position = [552 95 75 22];

            % Create SelectedMetabolitesPanel
            app.SelectedMetabolitesPanel = uipanel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.SelectedMetabolitesPanel.Tooltip = {'Change the metabolite list of the LCM'};
            app.SelectedMetabolitesPanel.ForegroundColor = [0.0392 0.2706 0.4314];
            app.SelectedMetabolitesPanel.BorderType = 'none';
            app.SelectedMetabolitesPanel.Title = 'Selected Metabolites';
            app.SelectedMetabolitesPanel.BackgroundColor = [1 0.9882 0.9882];
            app.SelectedMetabolitesPanel.FontWeight = 'bold';
            app.SelectedMetabolitesPanel.FontSize = 15;
            app.SelectedMetabolitesPanel.Position = [7 317 625 161];

            % Create AlaCheckBox
            app.AlaCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AlaCheckBox.Enable = 'off';
            app.AlaCheckBox.Text = 'Ala';
            app.AlaCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.AlaCheckBox.Position = [9 112 39 22];

            % Create AscCheckBox
            app.AscCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AscCheckBox.Enable = 'off';
            app.AscCheckBox.Text = 'Asc';
            app.AscCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.AscCheckBox.Position = [9 86 42 22];
            app.AscCheckBox.Value = true;

            % Create AspCheckBox
            app.AspCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AspCheckBox.Enable = 'off';
            app.AspCheckBox.Text = 'Asp';
            app.AspCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.AspCheckBox.Position = [9 60 43 22];
            app.AspCheckBox.Value = true;

            % Create bHBCheckBox
            app.bHBCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHBCheckBox.Enable = 'off';
            app.bHBCheckBox.Text = 'bHB';
            app.bHBCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.bHBCheckBox.Position = [9 34 46 22];

            % Create bHGCheckBox
            app.bHGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHGCheckBox.Enable = 'off';
            app.bHGCheckBox.Text = 'bHG';
            app.bHGCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.bHGCheckBox.Position = [9 8 47 22];

            % Create CitCheckBox
            app.CitCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CitCheckBox.Enable = 'off';
            app.CitCheckBox.Text = 'Cit';
            app.CitCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.CitCheckBox.Position = [66 112 37 22];

            % Create CrCheckBox
            app.CrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCheckBox.Enable = 'off';
            app.CrCheckBox.Text = 'Cr';
            app.CrCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.CrCheckBox.Position = [66 86 35 22];
            app.CrCheckBox.Value = true;

            % Create CrCH2CheckBox
            app.CrCH2CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCH2CheckBox.Enable = 'off';
            app.CrCH2CheckBox.Text = 'CrCH2';
            app.CrCH2CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.CrCH2CheckBox.Position = [67 60 59 22];
            app.CrCH2CheckBox.Value = true;

            % Create EtOHCheckBox
            app.EtOHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.EtOHCheckBox.Enable = 'off';
            app.EtOHCheckBox.Text = 'EtOH';
            app.EtOHCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.EtOHCheckBox.Position = [67 34 51 22];

            % Create GABACheckBox
            app.GABACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GABACheckBox.Enable = 'off';
            app.GABACheckBox.Text = 'GABA';
            app.GABACheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GABACheckBox.Position = [67 8 55 22];
            app.GABACheckBox.Value = true;

            % Create GPCCheckBox
            app.GPCCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GPCCheckBox.Enable = 'off';
            app.GPCCheckBox.Text = 'GPC';
            app.GPCCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GPCCheckBox.Position = [129 112 48 22];
            app.GPCCheckBox.Value = true;

            % Create GSHCheckBox
            app.GSHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GSHCheckBox.Enable = 'off';
            app.GSHCheckBox.Text = 'GSH';
            app.GSHCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GSHCheckBox.Position = [129 86 48 22];
            app.GSHCheckBox.Value = true;

            % Create GlcCheckBox
            app.GlcCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlcCheckBox.Enable = 'off';
            app.GlcCheckBox.Text = 'Glc';
            app.GlcCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GlcCheckBox.Position = [129 60 40 22];

            % Create GlnCheckBox
            app.GlnCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlnCheckBox.Enable = 'off';
            app.GlnCheckBox.Text = 'Gln';
            app.GlnCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GlnCheckBox.Position = [129 34 40 22];
            app.GlnCheckBox.Value = true;

            % Create GlyCheckBox
            app.GlyCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlyCheckBox.Enable = 'off';
            app.GlyCheckBox.Text = 'Gly';
            app.GlyCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GlyCheckBox.Position = [187 112 40 22];

            % Create H2OCheckBox
            app.H2OCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.H2OCheckBox.Enable = 'off';
            app.H2OCheckBox.Text = 'H2O';
            app.H2OCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.H2OCheckBox.Position = [187 85 46 22];
            app.H2OCheckBox.Value = true;

            % Create InsCheckBox
            app.InsCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.InsCheckBox.Enable = 'off';
            app.InsCheckBox.Text = 'Ins';
            app.InsCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.InsCheckBox.Position = [187 57 38 22];
            app.InsCheckBox.Value = true;

            % Create LacCheckBox
            app.LacCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.LacCheckBox.Enable = 'off';
            app.LacCheckBox.Text = 'Lac';
            app.LacCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.LacCheckBox.Position = [187 33 42 22];
            app.LacCheckBox.Value = true;

            % Create NAACheckBox
            app.NAACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAACheckBox.Enable = 'off';
            app.NAACheckBox.Text = 'NAA';
            app.NAACheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.NAACheckBox.Position = [187 7 46 22];
            app.NAACheckBox.Value = true;

            % Create GluCheckBox
            app.GluCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GluCheckBox.Enable = 'off';
            app.GluCheckBox.Text = 'Glu';
            app.GluCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.GluCheckBox.Position = [129 7 40 22];
            app.GluCheckBox.Value = true;

            % Create PChCheckBox
            app.PChCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PChCheckBox.Enable = 'off';
            app.PChCheckBox.Text = 'PCh';
            app.PChCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.PChCheckBox.Position = [250 85 45 22];
            app.PChCheckBox.Value = true;

            % Create PCrCheckBox
            app.PCrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PCrCheckBox.Enable = 'off';
            app.PCrCheckBox.Text = 'PCr';
            app.PCrCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.PCrCheckBox.Position = [250 59 42 22];
            app.PCrCheckBox.Value = true;

            % Create PECheckBox
            app.PECheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PECheckBox.Enable = 'off';
            app.PECheckBox.Text = 'PE';
            app.PECheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.PECheckBox.Position = [250 33 37 22];
            app.PECheckBox.Value = true;

            % Create PhenylCheckBox
            app.PhenylCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PhenylCheckBox.Enable = 'off';
            app.PhenylCheckBox.Text = 'Phenyl';
            app.PhenylCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.PhenylCheckBox.Position = [250 7 58 22];

            % Create ScylloCheckBox
            app.ScylloCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.ScylloCheckBox.Enable = 'off';
            app.ScylloCheckBox.Text = 'Scyllo';
            app.ScylloCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.ScylloCheckBox.Position = [320 112 54 22];
            app.ScylloCheckBox.Value = true;

            % Create SerCheckBox
            app.SerCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.SerCheckBox.Enable = 'off';
            app.SerCheckBox.Text = 'Ser';
            app.SerCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.SerCheckBox.Position = [320 86 40 22];

            % Create TauCheckBox
            app.TauCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TauCheckBox.Enable = 'off';
            app.TauCheckBox.Text = 'Tau';
            app.TauCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.TauCheckBox.Position = [320 60 41 22];
            app.TauCheckBox.Value = true;

            % Create NAAGCheckBox
            app.NAAGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAAGCheckBox.Enable = 'off';
            app.NAAGCheckBox.Text = 'NAAG';
            app.NAAGCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.NAAGCheckBox.Position = [251 112 55 22];
            app.NAAGCheckBox.Value = true;

            % Create TyrosCheckBox
            app.TyrosCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TyrosCheckBox.Enable = 'off';
            app.TyrosCheckBox.Text = 'Tyros';
            app.TyrosCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.TyrosCheckBox.Position = [320 33 50 22];

            % Create MM09CheckBox
            app.MM09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM09CheckBox.Enable = 'off';
            app.MM09CheckBox.Text = 'MM09';
            app.MM09CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM09CheckBox.Position = [388 110 56 22];
            app.MM09CheckBox.Value = true;

            % Create MM12CheckBox
            app.MM12CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM12CheckBox.Enable = 'off';
            app.MM12CheckBox.Text = 'MM12';
            app.MM12CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM12CheckBox.Position = [388 84 56 22];
            app.MM12CheckBox.Value = true;

            % Create MM14CheckBox
            app.MM14CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM14CheckBox.Enable = 'off';
            app.MM14CheckBox.Text = 'MM14';
            app.MM14CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM14CheckBox.Position = [388 58 56 22];
            app.MM14CheckBox.Value = true;

            % Create MM17CheckBox
            app.MM17CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM17CheckBox.Enable = 'off';
            app.MM17CheckBox.Text = 'MM17';
            app.MM17CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM17CheckBox.Position = [388 32 56 22];
            app.MM17CheckBox.Value = true;

            % Create MM20CheckBox
            app.MM20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM20CheckBox.Enable = 'off';
            app.MM20CheckBox.Text = 'MM20';
            app.MM20CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM20CheckBox.Position = [388 6 56 22];
            app.MM20CheckBox.Value = true;

            % Create Lip09CheckBox
            app.Lip09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip09CheckBox.Enable = 'off';
            app.Lip09CheckBox.Text = 'Lip09';
            app.Lip09CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.Lip09CheckBox.Position = [456 111 52 22];
            app.Lip09CheckBox.Value = true;

            % Create Lip13CheckBox
            app.Lip13CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip13CheckBox.Enable = 'off';
            app.Lip13CheckBox.Text = 'Lip13';
            app.Lip13CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.Lip13CheckBox.Position = [456 85 52 22];
            app.Lip13CheckBox.Value = true;

            % Create Lip20CheckBox
            app.Lip20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip20CheckBox.Enable = 'off';
            app.Lip20CheckBox.Text = 'Lip20';
            app.Lip20CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.Lip20CheckBox.Position = [456 59 52 22];
            app.Lip20CheckBox.Value = true;

            % Create MM37CheckBox
            app.MM37CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM37CheckBox.Enable = 'off';
            app.MM37CheckBox.Text = 'MM37';
            app.MM37CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM37CheckBox.Position = [456 32 56 22];

            % Create MM38CheckBox
            app.MM38CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM38CheckBox.Enable = 'off';
            app.MM38CheckBox.Text = 'MM38';
            app.MM38CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM38CheckBox.Position = [456 6 56 22];

            % Create MM40CheckBox
            app.MM40CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM40CheckBox.Enable = 'off';
            app.MM40CheckBox.Text = 'MM40';
            app.MM40CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM40CheckBox.Position = [523 111 56 22];

            % Create MM42CheckBox
            app.MM42CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM42CheckBox.Enable = 'off';
            app.MM42CheckBox.Text = 'MM42';
            app.MM42CheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM42CheckBox.Position = [522 86 56 22];

            % Create MMexpCheckBox
            app.MMexpCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MMexpCheckBox.Enable = 'off';
            app.MMexpCheckBox.Text = 'MMexp';
            app.MMexpCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MMexpCheckBox.Position = [522 60 63 22];

            % Create MM_PCCCheckBox
            app.MM_PCCCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM_PCCCheckBox.Enable = 'off';
            app.MM_PCCCheckBox.Text = 'MM_PCC';
            app.MM_PCCCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM_PCCCheckBox.Position = [522 33 74 22];

            % Create MM_CSOCheckBox
            app.MM_CSOCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM_CSOCheckBox.Enable = 'off';
            app.MM_CSOCheckBox.Text = 'MM_CSO';
            app.MM_CSOCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.MM_CSOCheckBox.Position = [522 7 74 22];

            % Create CystatCheckBox
            app.CystatCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CystatCheckBox.Enable = 'off';
            app.CystatCheckBox.Text = 'Cystat';
            app.CystatCheckBox.FontColor = [0.0392 0.2706 0.4314];
            app.CystatCheckBox.Position = [320 8 57 22];

            % Create SpecifyMRSandAnatomicalImagingFilesPanel
            app.SpecifyMRSandAnatomicalImagingFilesPanel = uipanel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.SpecifyMRSandAnatomicalImagingFilesPanel.Tooltip = {'Add MRS and anatomical image files here. Change the number of expected files first and use the SPM dialog to find the data.'};
            app.SpecifyMRSandAnatomicalImagingFilesPanel.ForegroundColor = [0.0392 0.2706 0.4314];
            app.SpecifyMRSandAnatomicalImagingFilesPanel.BorderType = 'none';
            app.SpecifyMRSandAnatomicalImagingFilesPanel.Title = '3. Specify MRS and Anatomical Imaging Files';
            app.SpecifyMRSandAnatomicalImagingFilesPanel.BackgroundColor = [1 0.9882 0.9882];
            app.SpecifyMRSandAnatomicalImagingFilesPanel.FontWeight = 'bold';
            app.SpecifyMRSandAnatomicalImagingFilesPanel.FontSize = 15;
            app.SpecifyMRSandAnatomicalImagingFilesPanel.Position = [6 144 627 169];

            % Create MRSDataButton
            app.MRSDataButton = uibutton(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'push');
            app.MRSDataButton.ButtonPushedFcn = createCallbackFcn(app, @MRSDataButtonPushed, true);
            app.MRSDataButton.BackgroundColor = [0.8 0.8 0.8];
            app.MRSDataButton.FontWeight = 'bold';
            app.MRSDataButton.FontColor = [0.0392 0.2706 0.4314];
            app.MRSDataButton.Tooltip = {'Add metaoblite data.'};
            app.MRSDataButton.Position = [10 72 118 23];
            app.MRSDataButton.Text = 'MRS Data';

            % Create H2OReferenceButton
            app.H2OReferenceButton = uibutton(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'push');
            app.H2OReferenceButton.ButtonPushedFcn = createCallbackFcn(app, @H2OReferenceButtonPushed, true);
            app.H2OReferenceButton.BackgroundColor = [0.8 0.8 0.8];
            app.H2OReferenceButton.FontWeight = 'bold';
            app.H2OReferenceButton.FontColor = [0.0392 0.2706 0.4314];
            app.H2OReferenceButton.Enable = 'off';
            app.H2OReferenceButton.Tooltip = {'Add eddy-current correction water reference (optional).'};
            app.H2OReferenceButton.Position = [321 72 118 23];
            app.H2OReferenceButton.Text = 'H2O Reference';

            % Create H2OShortTEButton
            app.H2OShortTEButton = uibutton(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'push');
            app.H2OShortTEButton.ButtonPushedFcn = createCallbackFcn(app, @H2OShortTEButtonPushed, true);
            app.H2OShortTEButton.BackgroundColor = [0.8 0.8 0.8];
            app.H2OShortTEButton.FontWeight = 'bold';
            app.H2OShortTEButton.FontColor = [0.0392 0.2706 0.4314];
            app.H2OShortTEButton.Enable = 'off';
            app.H2OShortTEButton.Tooltip = {'Add short-TE water reference for quantification (optional)'};
            app.H2OShortTEButton.Position = [9 42 118 23];
            app.H2OShortTEButton.Text = 'H2O Short TE';

            % Create MetaboliteNulledButton
            app.MetaboliteNulledButton = uibutton(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'push');
            app.MetaboliteNulledButton.ButtonPushedFcn = createCallbackFcn(app, @MetaboliteNulledButtonPushed, true);
            app.MetaboliteNulledButton.BackgroundColor = [0.8 0.8 0.8];
            app.MetaboliteNulledButton.FontWeight = 'bold';
            app.MetaboliteNulledButton.FontColor = [0.0392 0.2706 0.4314];
            app.MetaboliteNulledButton.Enable = 'off';
            app.MetaboliteNulledButton.Tooltip = {'Add metabolite-nulled spectra (optional).'};
            app.MetaboliteNulledButton.Position = [321 42 118 23];
            app.MetaboliteNulledButton.Text = 'Metabolite-Nulled';

            % Create T1DataniftiniiButton
            app.T1DataniftiniiButton = uibutton(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'push');
            app.T1DataniftiniiButton.ButtonPushedFcn = createCallbackFcn(app, @T1DataniftiniiButtonPushed, true);
            app.T1DataniftiniiButton.BackgroundColor = [0.8 0.8 0.8];
            app.T1DataniftiniiButton.FontWeight = 'bold';
            app.T1DataniftiniiButton.FontColor = [0.0392 0.2706 0.4314];
            app.T1DataniftiniiButton.Enable = 'off';
            app.T1DataniftiniiButton.Tooltip = {'Add T1-weighted anatomical scan for coregistration and segmentation (optional).'};
            app.T1DataniftiniiButton.Position = [10 6 117 23];
            app.T1DataniftiniiButton.Text = 'T1 Data (nifti *.nii)';

            % Create NumberofdatasetsEditFieldLabel
            app.NumberofdatasetsEditFieldLabel = uilabel(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.NumberofdatasetsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofdatasetsEditFieldLabel.FontColor = [0.0392 0.2706 0.4314];
            app.NumberofdatasetsEditFieldLabel.Position = [17 112 112 22];
            app.NumberofdatasetsEditFieldLabel.Text = 'Number of datasets';

            % Create NumberofdatasetsEditField
            app.NumberofdatasetsEditField = uieditfield(app.SpecifyMRSandAnatomicalImagingFilesPanel, 'numeric');
            app.NumberofdatasetsEditField.Limits = [1 Inf];
            app.NumberofdatasetsEditField.ValueDisplayFormat = '%.0f';
            app.NumberofdatasetsEditField.HorizontalAlignment = 'center';
            app.NumberofdatasetsEditField.FontColor = [0.0392 0.2706 0.4314];
            app.NumberofdatasetsEditField.Tooltip = {'Change the number of subjects.'};
            app.NumberofdatasetsEditField.Position = [144 112 100 22];
            app.NumberofdatasetsEditField.Value = 1;

            % Create MRSDataText
            app.MRSDataText = uitextarea(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.MRSDataText.Editable = 'off';
            app.MRSDataText.FontColor = [0.0392 0.2706 0.4314];
            app.MRSDataText.Position = [144 71 153 25];

            % Create H2OShortTEText
            app.H2OShortTEText = uitextarea(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.H2OShortTEText.Editable = 'off';
            app.H2OShortTEText.FontColor = [0.0392 0.2706 0.4314];
            app.H2OShortTEText.Position = [144 41 154 25];

            % Create H2OReferenceText
            app.H2OReferenceText = uitextarea(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.H2OReferenceText.Editable = 'off';
            app.H2OReferenceText.FontColor = [0.0392 0.2706 0.4314];
            app.H2OReferenceText.Position = [458 71 152 25];

            % Create MetaboliteNulledText
            app.MetaboliteNulledText = uitextarea(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.MetaboliteNulledText.Editable = 'off';
            app.MetaboliteNulledText.FontColor = [0.0392 0.2706 0.4314];
            app.MetaboliteNulledText.Position = [458 41 152 25];

            % Create T1DataText
            app.T1DataText = uitextarea(app.SpecifyMRSandAnatomicalImagingFilesPanel);
            app.T1DataText.Editable = 'off';
            app.T1DataText.FontColor = [0.0392 0.2706 0.4314];
            app.T1DataText.Position = [143 5 154 25];

            % Create SpecifyOutputFolderPanel
            app.SpecifyOutputFolderPanel = uipanel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.SpecifyOutputFolderPanel.Tooltip = {'Specify the output folder and the name of the jobfile.'};
            app.SpecifyOutputFolderPanel.ForegroundColor = [0.0392 0.2706 0.4314];
            app.SpecifyOutputFolderPanel.BorderType = 'none';
            app.SpecifyOutputFolderPanel.Title = '4. Specify Output Folder';
            app.SpecifyOutputFolderPanel.BackgroundColor = [1 0.9882 0.9882];
            app.SpecifyOutputFolderPanel.FontWeight = 'bold';
            app.SpecifyOutputFolderPanel.FontSize = 15;
            app.SpecifyOutputFolderPanel.Position = [6 68 627 72];

            % Create OutputFolderEditField
            app.OutputFolderEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.OutputFolderEditField.Editable = 'off';
            app.OutputFolderEditField.FontColor = [0.0392 0.2706 0.4314];
            app.OutputFolderEditField.Position = [144 14 183 22];

            % Create OutputFolderButton
            app.OutputFolderButton = uibutton(app.SpecifyOutputFolderPanel, 'push');
            app.OutputFolderButton.ButtonPushedFcn = createCallbackFcn(app, @OutputFolderButtonPushed, true);
            app.OutputFolderButton.BackgroundColor = [0.8 0.8 0.8];
            app.OutputFolderButton.FontWeight = 'bold';
            app.OutputFolderButton.FontColor = [0.0392 0.2706 0.4314];
            app.OutputFolderButton.Tooltip = {'Specifiy the output folder.'};
            app.OutputFolderButton.Position = [9 13 119 23];
            app.OutputFolderButton.Text = 'Output Folder';

            % Create JobNameEditFieldLabel
            app.JobNameEditFieldLabel = uilabel(app.SpecifyOutputFolderPanel);
            app.JobNameEditFieldLabel.HorizontalAlignment = 'right';
            app.JobNameEditFieldLabel.FontColor = [0.0392 0.2706 0.4314];
            app.JobNameEditFieldLabel.Position = [342 14 61 22];
            app.JobNameEditFieldLabel.Text = 'Job Name';

            % Create JobNameEditField
            app.JobNameEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.JobNameEditField.FontColor = [0.0392 0.2706 0.4314];
            app.JobNameEditField.Position = [468 14 142 22];
            app.JobNameEditField.Value = 'OspreyJob';

            % Create StatcsvEditField
            app.StatcsvEditField = uieditfield(app.InteractiveOspreyjobfilegeneratorUIFigure, 'text');
            app.StatcsvEditField.Editable = 'off';
            app.StatcsvEditField.FontColor = [0.0392 0.2706 0.4314];
            app.StatcsvEditField.Position = [150 47 183 22];

            % Create CANCELButton
            app.CANCELButton = uibutton(app.InteractiveOspreyjobfilegeneratorUIFigure, 'push');
            app.CANCELButton.ButtonPushedFcn = createCallbackFcn(app, @CANCELButtonPushed, true);
            app.CANCELButton.BackgroundColor = [1 0.9882 0.9882];
            app.CANCELButton.FontSize = 15;
            app.CANCELButton.FontWeight = 'bold';
            app.CANCELButton.FontColor = [0.0392 0.2706 0.4314];
            app.CANCELButton.Tooltip = {'See you next time!'};
            app.CANCELButton.Position = [314 8 318 32];
            app.CANCELButton.Text = 'CANCEL';

            % Create CREATEJOBButton
            app.CREATEJOBButton = uibutton(app.InteractiveOspreyjobfilegeneratorUIFigure, 'push');
            app.CREATEJOBButton.ButtonPushedFcn = createCallbackFcn(app, @CREATEJOBButtonPushed, true);
            app.CREATEJOBButton.BackgroundColor = [1 0.9882 0.9882];
            app.CREATEJOBButton.FontSize = 15;
            app.CREATEJOBButton.FontWeight = 'bold';
            app.CREATEJOBButton.FontColor = [0.0392 0.2706 0.4314];
            app.CREATEJOBButton.Tooltip = {'Generate jobfile.'};
            app.CREATEJOBButton.Position = [8 8 296 32];
            app.CREATEJOBButton.Text = 'CREATE JOB';

            % Create Image
            app.Image = uiimage(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.Image.Position = [590 822 67 29];
            app.Image.ImageSource = 'osprey.gif';

            % Create StatcsvFileButton
            app.StatcsvFileButton = uibutton(app.InteractiveOspreyjobfilegeneratorUIFigure, 'push');
            app.StatcsvFileButton.ButtonPushedFcn = createCallbackFcn(app, @StatcsvFileButtonPushed, true);
            app.StatcsvFileButton.BackgroundColor = [0.8 0.8 0.8];
            app.StatcsvFileButton.FontWeight = 'bold';
            app.StatcsvFileButton.FontColor = [0.0392 0.2706 0.4314];
            app.StatcsvFileButton.Tooltip = {'Specifiy the statistics csv file.'};
            app.StatcsvFileButton.Position = [15 46 119 23];
            app.StatcsvFileButton.Text = 'Stat csv File';

            % Create BasisSetEditField
            app.BasisSetEditField = uieditfield(app.InteractiveOspreyjobfilegeneratorUIFigure, 'text');
            app.BasisSetEditField.Editable = 'off';
            app.BasisSetEditField.FontColor = [0.0392 0.2706 0.4314];
            app.BasisSetEditField.Position = [453 482 183 22];

            % Create basissetfileButton
            app.basissetfileButton = uibutton(app.InteractiveOspreyjobfilegeneratorUIFigure, 'push');
            app.basissetfileButton.ButtonPushedFcn = createCallbackFcn(app, @basissetfileButtonPushed, true);
            app.basissetfileButton.BackgroundColor = [0.8 0.8 0.8];
            app.basissetfileButton.FontWeight = 'bold';
            app.basissetfileButton.FontColor = [0.0392 0.2706 0.4314];
            app.basissetfileButton.Tooltip = {'Specify basis set file in .mat format'};
            app.basissetfileButton.Position = [318 481 119 23];
            app.basissetfileButton.Text = 'basis set file';

            % Create FWHMMM3coLabel
            app.FWHMMM3coLabel = uilabel(app.InteractiveOspreyjobfilegeneratorUIFigure);
            app.FWHMMM3coLabel.HorizontalAlignment = 'right';
            app.FWHMMM3coLabel.FontColor = [0.0392 0.2706 0.4314];
            app.FWHMMM3coLabel.Enable = 'off';
            app.FWHMMM3coLabel.Position = [499 754 87 22];
            app.FWHMMM3coLabel.Text = 'FWHM MM3co';

            % Create FWHMMM3coEditField
            app.FWHMMM3coEditField = uieditfield(app.InteractiveOspreyjobfilegeneratorUIFigure, 'numeric');
            app.FWHMMM3coEditField.Limits = [0 Inf];
            app.FWHMMM3coEditField.HorizontalAlignment = 'left';
            app.FWHMMM3coEditField.FontColor = [0.0392 0.2706 0.4314];
            app.FWHMMM3coEditField.Enable = 'off';
            app.FWHMMM3coEditField.Tooltip = {'FWHM of the co-edited MM peak at 3 ppm.'};
            app.FWHMMM3coEditField.Position = [601 754 33 22];
            app.FWHMMM3coEditField.Value = 14;

            % Show the figure after all components are created
            app.InteractiveOspreyjobfilegeneratorUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CreateOspreyJob_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.InteractiveOspreyjobfilegeneratorUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.InteractiveOspreyjobfilegeneratorUIFigure)
        end
    end
end