classdef CreateOspreyJob_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        InteractiveOspreyJobmCreatorUIFigure  matlab.ui.Figure
        SpecifySequenceInformationPanel  matlab.ui.container.Panel
        SequenceTypeDropDownLabel       matlab.ui.control.Label
        SequenceTypeDropDown            matlab.ui.control.DropDown
        EditingTargetsDropDownLabel     matlab.ui.control.Label
        EditingTargetsDropDown          matlab.ui.control.DropDown
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
        MRSProtocolDropDownLabel        matlab.ui.control.Label
        MRSProtocolDropDown             matlab.ui.control.DropDown
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
        MRSdataDropDown                 matlab.ui.control.DropDown
        H2OReferenceDropDown            matlab.ui.control.DropDown
        H2OShortTEDropDown              matlab.ui.control.DropDown
        MetaboliteNulledDropDown        matlab.ui.control.DropDown
        T1DataniftiniiDropDown          matlab.ui.control.DropDown
        SpecifyOutputFolderPanel        matlab.ui.container.Panel
        OutputFolderEditField           matlab.ui.control.EditField
        OutputFolderButton              matlab.ui.control.Button
        JobNameEditFieldLabel           matlab.ui.control.Label
        JobNameEditField                matlab.ui.control.EditField
        ExitButton                      matlab.ui.control.Button
        RunJobButton                    matlab.ui.control.Button
        LoadJobFileButton               matlab.ui.control.Button
        LoadMRSContButton               matlab.ui.control.Button
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
                    app.EditingTargetsDropDown.Items = {'GABA','GSH'};
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
            info = 'Please select the MRS file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            mrsfiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            if isempty(mrsfiles)
                return
            end

            [~,file_basename,file_exten]=fileparts(mrsfiles(1,:));
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} mrsfiles(i,:)};
            end
            
            app.MRSdataDropDown.Items = filelist;
            app.MRSdataDropDown.Value = filelist{1};
            
            if file_exten=='.7'
                app.H2OReferenceButton.Enable = 'Off';
                app.H2OReferenceDropDown.Enable = 'Off';
                app.H2OReferenceDropDown.Items = {};
                app.H2OReferenceDropDown.Value = {};
            else
                app.H2OReferenceButton.Enable = 'On';
                app.H2OReferenceDropDown.Enable = 'On';
                app.H2OReferenceDropDown.Items = {};
                app.H2OReferenceDropDown.Value = {};
            end
            
            app.H2OShortTEButton.Enable = 'On';
            app.H2OShortTEDropDown.Enable = 'On';
            app.MetaboliteNulledButton.Enable = 'On';
            app.MetaboliteNulledDropDown.Enable = 'On';
            app.T1DataniftiniiButton.Enable = 'On';
            app.T1DataniftiniiDropDown.Enable = 'On';
        end

        % Button pushed function: H2OReferenceButton
        function H2OReferenceButtonPushed(app, event)
            info = 'Please select the water reference file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            h2oreffiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            if isempty(h2oreffiles)
                return
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2oreffiles(i,:)};
            end
            
            app.H2OReferenceDropDown.Items = filelist;
            app.H2OReferenceDropDown.Value = filelist{1};
            
        end

        % Button pushed function: H2OShortTEButton
        function H2OShortTEButtonPushed(app, event)
            info = 'Please select the water short-TE file to read';
            [fname,pathname]=uigetfile('*.*',info);
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            h2ostefiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            if isempty(h2ostefiles)
                return
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} h2ostefiles(i,:)};
            end
            
            app.H2OShortTEDropDown.Items = filelist;
            app.H2OShortTEDropDown.Value = filelist{1};
        end

        % Button pushed function: MetaboliteNulledButton
        function MetaboliteNulledButtonPushed(app, event)
            info = 'Please select the metabolite-nulled file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            metnulfiles = spm_select(ndata,'any',info,{},pwd,'.*','1');
            
            if isempty(metnulfiles)
                return
            end
            
            filelist = {};
            for i=1:ndata
                filelist = {filelist{:} metnulfiles(i,:)};
            end
            
            app.MetaboliteNulledDropDown.Items = filelist;
            app.MetaboliteNulledDropDown.Value = filelist{1};
        end

        % Button pushed function: T1DataniftiniiButton
        function T1DataniftiniiButtonPushed(app, event)
            info = 'Please select the T1 Data file to read';
            
            ndata = app.NumberofdatasetsEditField.Value;
            
            t1imfiles = spm_select(ndata,'image',info,{},pwd,'.nii','1');
            
            if isempty(t1imfiles)
                return
            end
            
            filelist = {};
            for i=1:ndata
                splniftiFile = split(t1imfiles(i,:),',');
                filelist = {filelist{:} splniftiFile{1}};
            end
            
            app.T1DataniftiniiDropDown.Items = filelist;
            app.T1DataniftiniiDropDown.Value = filelist{1};
        end

        % Button pushed function: OutputFolderButton
        function OutputFolderButtonPushed(app, event)
            info = 'Please select the output folder';
            pathname=uigetdir('*.*',info);
            
            if ~pathname
                return
            end
            
            app.OutputFolderEditField.Value = pathname;
        end

        % Button pushed function: ExitButton
        function ExitButtonPushed(app, event)
            [ospFFolder,~,~] = fileparts(which('Osprey.m'));
            curdir = cd(ospFFolder);
            load('startpath.mat')
            path(startpath);
            delete('startpath.mat')
            cd(curdir)
            
            delete(app.InteractiveOspreyJobmCreatorUIFigure)
        end

        % Button pushed function: RunJobButton
        function RunJobButtonPushed(app, event)
            jobm = osp_create_job_file(app);
            
            delete(app.InteractiveOspreyJobmCreatorUIFigure)
            
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

        % Button pushed function: LoadJobFileButton
        function LoadJobFileButtonPushed(app, event)
            info = 'Select Osprey job file (.m or .csv-file)';
            [mfile,mpath] = uigetfile({'*.m','*.csv'},info);
            
            if ~mfile
                return
            end
            
            [~,~,file_exten]=fileparts(fullfile(mpath,mfile));
            
            switch file_exten
                case '.m'
                    osp_set_osppreyapp_parameters(mfile,mpath,app);
                case '.csv'
                    delete(app.InteractiveOspreyJobmCreatorUIFigure);
                    
                    gui.data.MRSCont = OspreyJob(fullfile(mpath,mfile),1);
                    
                    %Default colormap
                    gui.colormap.Background = [255/255 254/255 254/255];
                    gui.colormap.LightAccent = [110/255 136/255 164/255];
                    gui.colormap.Foreground = [11/255 71/255 111/255];
                    gui.colormap.Accent = [254/255 186/255 47/255];
                    
                    gui.data.MRSCont.colormap = gui.colormap;
                    gui.data.MRSCont.colormapidx = 1;
                    gui.data.MRSCont.flags.isToolChecked = 1;
                    OspreyGUI(gui.data.MRSCont);
            end
            
        end

        % Button pushed function: LoadMRSContButton
        function LoadMRSContButtonPushed(app, event)
            info = 'Select Osprey MRSCont file (.mat-file)';
            [mfile,mpath] = uigetfile('*.mat',info);
            
            if ~mfile
                return
            end
            
            delete(app.InteractiveOspreyJobmCreatorUIFigure);
            
            gui.data = load(fullfile(mpath,mfile), 'MRSCont');
            
            %Default colormap
            gui.colormap.Background = [255/255 254/255 254/255];
            gui.colormap.LightAccent = [110/255 136/255 164/255];
            gui.colormap.Foreground = [11/255 71/255 111/255];
            gui.colormap.Accent = [254/255 186/255 47/255];
            
            gui.data.MRSCont.colormap = gui.colormap;
                    gui.data.MRSCont.colormapidx = 1;
            gui.data.MRSCont.flags.isToolChecked = 1;
            OspreyGUI(gui.data.MRSCont);
        end

        % Value changed function: NumberofdatasetsEditField
        function NumberofdatasetsEditFieldValueChanged(app, event)
            value = app.NumberofdatasetsEditField.Value;
            
            app.MRSdataDropDown.Items = {};
            app.MRSdataDropDown.Value = {};
            
            app.H2OReferenceDropDown.Items = {};
            app.H2OReferenceDropDown.Value = {};
            app.H2OReferenceDropDown.Enable = 'Off';
            app.H2OReferenceButton.Enable = 'Off';
            
            app.H2OShortTEDropDown.Items = {};
            app.H2OShortTEDropDown.Value = {};
            app.H2OShortTEDropDown.Enable = 'Off';
            app.H2OShortTEButton.Enable = 'Off';
            
            app.MetaboliteNulledDropDown.Items = {};
            app.MetaboliteNulledDropDown.Value = {};
            app.MetaboliteNulledDropDown.Enable = 'Off';
            app.MetaboliteNulledButton.Enable = 'Off';
            
            app.T1DataniftiniiDropDown.Items = {};
            app.T1DataniftiniiDropDown.Value = {};
            app.T1DataniftiniiDropDown.Enable = 'Off';
            app.T1DataniftiniiButton.Enable = 'Off';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create InteractiveOspreyJobmCreatorUIFigure and hide until all components are created
            app.InteractiveOspreyJobmCreatorUIFigure = uifigure('Visible', 'off');
            app.InteractiveOspreyJobmCreatorUIFigure.Color = [1 1 1];
            app.InteractiveOspreyJobmCreatorUIFigure.Colormap = [0.2431 0.149 0.6588;0.2431 0.1529 0.6745;0.2471 0.1569 0.6863;0.2471 0.1608 0.698;0.251 0.1647 0.7059;0.251 0.1686 0.7176;0.2549 0.1725 0.7294;0.2549 0.1765 0.7412;0.2588 0.1804 0.749;0.2588 0.1843 0.7608;0.2627 0.1882 0.7725;0.2627 0.1922 0.7843;0.2627 0.1961 0.7922;0.2667 0.2 0.8039;0.2667 0.2039 0.8157;0.2706 0.2078 0.8235;0.2706 0.2157 0.8353;0.2706 0.2196 0.8431;0.2745 0.2235 0.851;0.2745 0.2275 0.8627;0.2745 0.2314 0.8706;0.2745 0.2392 0.8784;0.2784 0.2431 0.8824;0.2784 0.2471 0.8902;0.2784 0.2549 0.898;0.2784 0.2588 0.902;0.2784 0.2667 0.9098;0.2784 0.2706 0.9137;0.2784 0.2745 0.9216;0.2824 0.2824 0.9255;0.2824 0.2863 0.9294;0.2824 0.2941 0.9333;0.2824 0.298 0.9412;0.2824 0.3059 0.9451;0.2824 0.3098 0.949;0.2824 0.3137 0.9529;0.2824 0.3216 0.9569;0.2824 0.3255 0.9608;0.2824 0.3294 0.9647;0.2784 0.3373 0.9686;0.2784 0.3412 0.9686;0.2784 0.349 0.9725;0.2784 0.3529 0.9765;0.2784 0.3569 0.9804;0.2784 0.3647 0.9804;0.2745 0.3686 0.9843;0.2745 0.3765 0.9843;0.2745 0.3804 0.9882;0.2706 0.3843 0.9882;0.2706 0.3922 0.9922;0.2667 0.3961 0.9922;0.2627 0.4039 0.9922;0.2627 0.4078 0.9961;0.2588 0.4157 0.9961;0.2549 0.4196 0.9961;0.251 0.4275 0.9961;0.2471 0.4314 1;0.2431 0.4392 1;0.2353 0.4431 1;0.2314 0.451 1;0.2235 0.4549 1;0.2196 0.4627 0.9961;0.2118 0.4667 0.9961;0.2078 0.4745 0.9922;0.2 0.4784 0.9922;0.1961 0.4863 0.9882;0.1922 0.4902 0.9882;0.1882 0.498 0.9843;0.1843 0.502 0.9804;0.1843 0.5098 0.9804;0.1804 0.5137 0.9765;0.1804 0.5176 0.9725;0.1804 0.5255 0.9725;0.1804 0.5294 0.9686;0.1765 0.5333 0.9647;0.1765 0.5412 0.9608;0.1765 0.5451 0.9569;0.1765 0.549 0.9529;0.1765 0.5569 0.949;0.1725 0.5608 0.9451;0.1725 0.5647 0.9412;0.1686 0.5686 0.9373;0.1647 0.5765 0.9333;0.1608 0.5804 0.9294;0.1569 0.5843 0.9255;0.1529 0.5922 0.9216;0.1529 0.5961 0.9176;0.149 0.6 0.9137;0.149 0.6039 0.9098;0.1451 0.6078 0.9098;0.1451 0.6118 0.9059;0.1412 0.6196 0.902;0.1412 0.6235 0.898;0.1373 0.6275 0.898;0.1373 0.6314 0.8941;0.1333 0.6353 0.8941;0.1294 0.6392 0.8902;0.1255 0.6471 0.8902;0.1216 0.651 0.8863;0.1176 0.6549 0.8824;0.1137 0.6588 0.8824;0.1137 0.6627 0.8784;0.1098 0.6667 0.8745;0.1059 0.6706 0.8706;0.102 0.6745 0.8667;0.098 0.6784 0.8627;0.0902 0.6824 0.8549;0.0863 0.6863 0.851;0.0784 0.6902 0.8471;0.0706 0.6941 0.8392;0.0627 0.698 0.8353;0.0549 0.702 0.8314;0.0431 0.702 0.8235;0.0314 0.7059 0.8196;0.0235 0.7098 0.8118;0.0157 0.7137 0.8078;0.0078 0.7176 0.8;0.0039 0.7176 0.7922;0 0.7216 0.7882;0 0.7255 0.7804;0 0.7294 0.7765;0.0039 0.7294 0.7686;0.0078 0.7333 0.7608;0.0157 0.7333 0.7569;0.0235 0.7373 0.749;0.0353 0.7412 0.7412;0.051 0.7412 0.7373;0.0627 0.7451 0.7294;0.0784 0.7451 0.7216;0.0902 0.749 0.7137;0.102 0.7529 0.7098;0.1137 0.7529 0.702;0.1255 0.7569 0.6941;0.1373 0.7569 0.6863;0.1451 0.7608 0.6824;0.1529 0.7608 0.6745;0.1608 0.7647 0.6667;0.1686 0.7647 0.6588;0.1725 0.7686 0.651;0.1804 0.7686 0.6471;0.1843 0.7725 0.6392;0.1922 0.7725 0.6314;0.1961 0.7765 0.6235;0.2 0.7804 0.6157;0.2078 0.7804 0.6078;0.2118 0.7843 0.6;0.2196 0.7843 0.5882;0.2235 0.7882 0.5804;0.2314 0.7882 0.5725;0.2392 0.7922 0.5647;0.251 0.7922 0.5529;0.2588 0.7922 0.5451;0.2706 0.7961 0.5373;0.2824 0.7961 0.5255;0.2941 0.7961 0.5176;0.3059 0.8 0.5059;0.3176 0.8 0.498;0.3294 0.8 0.4863;0.3412 0.8 0.4784;0.3529 0.8 0.4667;0.3686 0.8039 0.4549;0.3804 0.8039 0.4471;0.3922 0.8039 0.4353;0.4039 0.8039 0.4235;0.4196 0.8039 0.4118;0.4314 0.8039 0.4;0.4471 0.8039 0.3922;0.4627 0.8 0.3804;0.4745 0.8 0.3686;0.4902 0.8 0.3569;0.5059 0.8 0.349;0.5176 0.8 0.3373;0.5333 0.7961 0.3255;0.5451 0.7961 0.3176;0.5608 0.7961 0.3059;0.5765 0.7922 0.2941;0.5882 0.7922 0.2824;0.6039 0.7882 0.2745;0.6157 0.7882 0.2627;0.6314 0.7843 0.251;0.6431 0.7843 0.2431;0.6549 0.7804 0.2314;0.6706 0.7804 0.2235;0.6824 0.7765 0.2157;0.698 0.7765 0.2078;0.7098 0.7725 0.2;0.7216 0.7686 0.1922;0.7333 0.7686 0.1843;0.7451 0.7647 0.1765;0.7608 0.7647 0.1725;0.7725 0.7608 0.1647;0.7843 0.7569 0.1608;0.7961 0.7569 0.1569;0.8078 0.7529 0.1529;0.8157 0.749 0.1529;0.8275 0.749 0.1529;0.8392 0.7451 0.1529;0.851 0.7451 0.1569;0.8588 0.7412 0.1569;0.8706 0.7373 0.1608;0.8824 0.7373 0.1647;0.8902 0.7373 0.1686;0.902 0.7333 0.1765;0.9098 0.7333 0.1804;0.9176 0.7294 0.1882;0.9255 0.7294 0.1961;0.9373 0.7294 0.2078;0.9451 0.7294 0.2157;0.9529 0.7294 0.2235;0.9608 0.7294 0.2314;0.9686 0.7294 0.2392;0.9765 0.7294 0.2431;0.9843 0.7333 0.2431;0.9882 0.7373 0.2431;0.9961 0.7412 0.2392;0.9961 0.7451 0.2353;0.9961 0.7529 0.2314;0.9961 0.7569 0.2275;0.9961 0.7608 0.2235;0.9961 0.7686 0.2196;0.9961 0.7725 0.2157;0.9961 0.7804 0.2078;0.9961 0.7843 0.2039;0.9961 0.7922 0.2;0.9922 0.7961 0.1961;0.9922 0.8039 0.1922;0.9922 0.8078 0.1922;0.9882 0.8157 0.1882;0.9843 0.8235 0.1843;0.9843 0.8275 0.1804;0.9804 0.8353 0.1804;0.9765 0.8392 0.1765;0.9765 0.8471 0.1725;0.9725 0.851 0.1686;0.9686 0.8588 0.1647;0.9686 0.8667 0.1647;0.9647 0.8706 0.1608;0.9647 0.8784 0.1569;0.9608 0.8824 0.1569;0.9608 0.8902 0.1529;0.9608 0.898 0.149;0.9608 0.902 0.149;0.9608 0.9098 0.1451;0.9608 0.9137 0.1412;0.9608 0.9216 0.1373;0.9608 0.9255 0.1333;0.9608 0.9333 0.1294;0.9647 0.9373 0.1255;0.9647 0.9451 0.1216;0.9647 0.949 0.1176;0.9686 0.9569 0.1098;0.9686 0.9608 0.1059;0.9725 0.9686 0.102;0.9725 0.9725 0.0941;0.9765 0.9765 0.0863;0.9765 0.9843 0.0824;1 1 1];
            app.InteractiveOspreyJobmCreatorUIFigure.Position = [100 100 722 636];
            app.InteractiveOspreyJobmCreatorUIFigure.Name = 'Interactive OspreyJob.m Creator';

            % Create SpecifySequenceInformationPanel
            app.SpecifySequenceInformationPanel = uipanel(app.InteractiveOspreyJobmCreatorUIFigure);
            app.SpecifySequenceInformationPanel.ForegroundColor = [0.0392 0.2784 0.4392];
            app.SpecifySequenceInformationPanel.Title = '1. Specify Sequence Information';
            app.SpecifySequenceInformationPanel.BackgroundColor = [1 1 1];
            app.SpecifySequenceInformationPanel.FontWeight = 'bold';
            app.SpecifySequenceInformationPanel.FontSize = 15;
            app.SpecifySequenceInformationPanel.Position = [8 524 350 102];

            % Create SequenceTypeDropDownLabel
            app.SequenceTypeDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.SequenceTypeDropDownLabel.BackgroundColor = [1 1 1];
            app.SequenceTypeDropDownLabel.HorizontalAlignment = 'right';
            app.SequenceTypeDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.SequenceTypeDropDownLabel.Position = [28 45 88 22];
            app.SequenceTypeDropDownLabel.Text = 'Sequence Type';

            % Create SequenceTypeDropDown
            app.SequenceTypeDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.SequenceTypeDropDown.Items = {'unedited', 'MEGA', 'HERMES', 'HERCULES'};
            app.SequenceTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SequenceTypeDropDownValueChanged, true);
            app.SequenceTypeDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.SequenceTypeDropDown.BackgroundColor = [1 1 1];
            app.SequenceTypeDropDown.Position = [131 45 129 22];
            app.SequenceTypeDropDown.Value = 'unedited';

            % Create EditingTargetsDropDownLabel
            app.EditingTargetsDropDownLabel = uilabel(app.SpecifySequenceInformationPanel);
            app.EditingTargetsDropDownLabel.BackgroundColor = [1 1 1];
            app.EditingTargetsDropDownLabel.HorizontalAlignment = 'right';
            app.EditingTargetsDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.EditingTargetsDropDownLabel.Position = [28 12 85 22];
            app.EditingTargetsDropDownLabel.Text = 'Editing Targets';

            % Create EditingTargetsDropDown
            app.EditingTargetsDropDown = uidropdown(app.SpecifySequenceInformationPanel);
            app.EditingTargetsDropDown.Items = {'none'};
            app.EditingTargetsDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.EditingTargetsDropDown.BackgroundColor = [1 1 1];
            app.EditingTargetsDropDown.Position = [131 12 129 22];
            app.EditingTargetsDropDown.Value = 'none';

            % Create SpecifyDataHandlingandModelingOptionsPanel
            app.SpecifyDataHandlingandModelingOptionsPanel = uipanel(app.InteractiveOspreyJobmCreatorUIFigure);
            app.SpecifyDataHandlingandModelingOptionsPanel.ForegroundColor = [0.0392 0.2784 0.4392];
            app.SpecifyDataHandlingandModelingOptionsPanel.Title = '2. Specify Data Handling and Modeling Options';
            app.SpecifyDataHandlingandModelingOptionsPanel.BackgroundColor = [1 1 1];
            app.SpecifyDataHandlingandModelingOptionsPanel.FontWeight = 'bold';
            app.SpecifyDataHandlingandModelingOptionsPanel.FontSize = 15;
            app.SpecifyDataHandlingandModelingOptionsPanel.Position = [8 167 350 343];

            % Create SaveLCMCheckBox
            app.SaveLCMCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveLCMCheckBox.Text = 'Save LCM';
            app.SaveLCMCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.SaveLCMCheckBox.Position = [21 283 78 22];

            % Create SaveJMRUICheckBox
            app.SaveJMRUICheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveJMRUICheckBox.Text = 'Save JMRUI';
            app.SaveJMRUICheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.SaveJMRUICheckBox.Position = [124 283 89 22];

            % Create SaveVendorCheckBox
            app.SaveVendorCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.SaveVendorCheckBox.Text = 'Save Vendor';
            app.SaveVendorCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.SaveVendorCheckBox.Position = [238 283 90 22];

            % Create FittingAlgorithmDropDownLabel
            app.FittingAlgorithmDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingAlgorithmDropDownLabel.BackgroundColor = [1 1 1];
            app.FittingAlgorithmDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.FittingAlgorithmDropDownLabel.Position = [28 213 94 22];
            app.FittingAlgorithmDropDownLabel.Text = 'Fitting Algorithm';

            % Create FittingAlgorithmDropDown
            app.FittingAlgorithmDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingAlgorithmDropDown.Items = {'Osprey', 'AQSES', 'TARQUIN'};
            app.FittingAlgorithmDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.FittingAlgorithmDropDown.BackgroundColor = [1 1 1];
            app.FittingAlgorithmDropDown.Position = [161 213 129 22];
            app.FittingAlgorithmDropDown.Value = 'Osprey';

            % Create IncludedMetabolitesDropDownLabel
            app.IncludedMetabolitesDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.IncludedMetabolitesDropDownLabel.BackgroundColor = [1 1 1];
            app.IncludedMetabolitesDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.IncludedMetabolitesDropDownLabel.Position = [28 180 118 22];
            app.IncludedMetabolitesDropDownLabel.Text = 'Included Metabolites';

            % Create IncludedMetabolitesDropDown
            app.IncludedMetabolitesDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.IncludedMetabolitesDropDown.Items = {'default', 'full', 'custom'};
            app.IncludedMetabolitesDropDown.ValueChangedFcn = createCallbackFcn(app, @IncludedMetabolitesDropDownValueChanged, true);
            app.IncludedMetabolitesDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.IncludedMetabolitesDropDown.BackgroundColor = [1 1 1];
            app.IncludedMetabolitesDropDown.Position = [161 180 129 22];
            app.IncludedMetabolitesDropDown.Value = 'default';

            % Create FittingStyleDropDownLabel
            app.FittingStyleDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingStyleDropDownLabel.BackgroundColor = [1 1 1];
            app.FittingStyleDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.FittingStyleDropDownLabel.Enable = 'off';
            app.FittingStyleDropDownLabel.Position = [28 146 69 22];
            app.FittingStyleDropDownLabel.Text = 'Fitting Style';

            % Create FittingStyleDropDown
            app.FittingStyleDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.FittingStyleDropDown.Items = {'Concatenated', 'Separate'};
            app.FittingStyleDropDown.Enable = 'off';
            app.FittingStyleDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.FittingStyleDropDown.BackgroundColor = [1 1 1];
            app.FittingStyleDropDown.Position = [161 146 129 22];
            app.FittingStyleDropDown.Value = 'Concatenated';

            % Create MRSFitRangeppmEditFieldLabel
            app.MRSFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.MRSFitRangeppmEditFieldLabel.FontColor = [0.0392 0.2784 0.4392];
            app.MRSFitRangeppmEditFieldLabel.Position = [28 113 121 22];
            app.MRSFitRangeppmEditFieldLabel.Text = 'MRS Fit Range (ppm)';

            % Create MRSFitRangeppmEditField
            app.MRSFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.MRSFitRangeppmEditField.FontColor = [0.0392 0.2784 0.4392];
            app.MRSFitRangeppmEditField.Position = [161 113 129 22];
            app.MRSFitRangeppmEditField.Value = '0.2 4.2';

            % Create WaterFitRangeppmEditFieldLabel
            app.WaterFitRangeppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.WaterFitRangeppmEditFieldLabel.FontColor = [0.0392 0.2784 0.4392];
            app.WaterFitRangeppmEditFieldLabel.Position = [28 80 126 22];
            app.WaterFitRangeppmEditFieldLabel.Text = 'Water Fit Range (ppm)';

            % Create WaterFitRangeppmEditField
            app.WaterFitRangeppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'text');
            app.WaterFitRangeppmEditField.FontColor = [0.0392 0.2784 0.4392];
            app.WaterFitRangeppmEditField.Position = [161 80 129 22];
            app.WaterFitRangeppmEditField.Value = '2.0 7.4';

            % Create BaselineknotspacingppmEditFieldLabel
            app.BaselineknotspacingppmEditFieldLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.BaselineknotspacingppmEditFieldLabel.FontColor = [0.0392 0.2784 0.4392];
            app.BaselineknotspacingppmEditFieldLabel.Position = [28 47 158 22];
            app.BaselineknotspacingppmEditFieldLabel.Text = 'Baseline knot spacing (ppm)';

            % Create BaselineknotspacingppmEditField
            app.BaselineknotspacingppmEditField = uieditfield(app.SpecifyDataHandlingandModelingOptionsPanel, 'numeric');
            app.BaselineknotspacingppmEditField.ValueDisplayFormat = '%11.1g';
            app.BaselineknotspacingppmEditField.HorizontalAlignment = 'left';
            app.BaselineknotspacingppmEditField.FontColor = [0.0392 0.2784 0.4392];
            app.BaselineknotspacingppmEditField.Position = [201 47 89 22];
            app.BaselineknotspacingppmEditField.Value = 0.4;

            % Create AddMMandLipbasisfunctionstofitCheckBox
            app.AddMMandLipbasisfunctionstofitCheckBox = uicheckbox(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.AddMMandLipbasisfunctionstofitCheckBox.Text = 'Add MM and Lip basis functions to fit';
            app.AddMMandLipbasisfunctionstofitCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.AddMMandLipbasisfunctionstofitCheckBox.Position = [28 15 225 22];
            app.AddMMandLipbasisfunctionstofitCheckBox.Value = true;

            % Create MRSProtocolDropDownLabel
            app.MRSProtocolDropDownLabel = uilabel(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.MRSProtocolDropDownLabel.BackgroundColor = [1 1 1];
            app.MRSProtocolDropDownLabel.FontColor = [0.0392 0.2784 0.4392];
            app.MRSProtocolDropDownLabel.Position = [28 248 85 22];
            app.MRSProtocolDropDownLabel.Text = 'MRS Protocol';

            % Create MRSProtocolDropDown
            app.MRSProtocolDropDown = uidropdown(app.SpecifyDataHandlingandModelingOptionsPanel);
            app.MRSProtocolDropDown.Items = {'Brain', 'Phantom'};
            app.MRSProtocolDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.MRSProtocolDropDown.BackgroundColor = [1 1 1];
            app.MRSProtocolDropDown.Position = [161 248 129 22];
            app.MRSProtocolDropDown.Value = 'Brain';

            % Create SelectedMetabolitesPanel
            app.SelectedMetabolitesPanel = uipanel(app.InteractiveOspreyJobmCreatorUIFigure);
            app.SelectedMetabolitesPanel.ForegroundColor = [0.0392 0.2784 0.4392];
            app.SelectedMetabolitesPanel.TitlePosition = 'centertop';
            app.SelectedMetabolitesPanel.Title = 'Selected Metabolites';
            app.SelectedMetabolitesPanel.BackgroundColor = [1 1 1];
            app.SelectedMetabolitesPanel.FontWeight = 'bold';
            app.SelectedMetabolitesPanel.FontSize = 15;
            app.SelectedMetabolitesPanel.Position = [8 9 707 146];

            % Create AlaCheckBox
            app.AlaCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AlaCheckBox.Enable = 'off';
            app.AlaCheckBox.Text = 'Ala';
            app.AlaCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.AlaCheckBox.Position = [15 89 39 22];

            % Create AscCheckBox
            app.AscCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AscCheckBox.Enable = 'off';
            app.AscCheckBox.Text = 'Asc';
            app.AscCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.AscCheckBox.Position = [75 89 42 22];
            app.AscCheckBox.Value = true;

            % Create AspCheckBox
            app.AspCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.AspCheckBox.Enable = 'off';
            app.AspCheckBox.Text = 'Asp';
            app.AspCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.AspCheckBox.Position = [134 89 43 22];
            app.AspCheckBox.Value = true;

            % Create bHBCheckBox
            app.bHBCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHBCheckBox.Enable = 'off';
            app.bHBCheckBox.Text = 'bHB';
            app.bHBCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.bHBCheckBox.Position = [194 89 46 22];

            % Create bHGCheckBox
            app.bHGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.bHGCheckBox.Enable = 'off';
            app.bHGCheckBox.Text = 'bHG';
            app.bHGCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.bHGCheckBox.Position = [255 89 47 22];

            % Create CitCheckBox
            app.CitCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CitCheckBox.Enable = 'off';
            app.CitCheckBox.Text = 'Cit';
            app.CitCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.CitCheckBox.Position = [316 89 37 22];

            % Create CrCheckBox
            app.CrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCheckBox.Enable = 'off';
            app.CrCheckBox.Text = 'Cr';
            app.CrCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.CrCheckBox.Position = [373 89 35 22];
            app.CrCheckBox.Value = true;

            % Create CrCH2CheckBox
            app.CrCH2CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.CrCH2CheckBox.Enable = 'off';
            app.CrCH2CheckBox.Text = 'CrCH2';
            app.CrCH2CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.CrCH2CheckBox.Position = [427 89 59 22];
            app.CrCH2CheckBox.Value = true;

            % Create EtOHCheckBox
            app.EtOHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.EtOHCheckBox.Enable = 'off';
            app.EtOHCheckBox.Text = 'EtOH';
            app.EtOHCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.EtOHCheckBox.Position = [486 89 51 22];

            % Create GABACheckBox
            app.GABACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GABACheckBox.Enable = 'off';
            app.GABACheckBox.Text = 'GABA';
            app.GABACheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GABACheckBox.Position = [545 89 55 22];
            app.GABACheckBox.Value = true;

            % Create GPCCheckBox
            app.GPCCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GPCCheckBox.Enable = 'off';
            app.GPCCheckBox.Text = 'GPC';
            app.GPCCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GPCCheckBox.Position = [604 89 48 22];
            app.GPCCheckBox.Value = true;

            % Create GSHCheckBox
            app.GSHCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GSHCheckBox.Enable = 'off';
            app.GSHCheckBox.Text = 'GSH';
            app.GSHCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GSHCheckBox.Position = [656 89 48 22];
            app.GSHCheckBox.Value = true;

            % Create GlcCheckBox
            app.GlcCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlcCheckBox.Enable = 'off';
            app.GlcCheckBox.Text = 'Glc';
            app.GlcCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GlcCheckBox.Position = [15 62 40 22];

            % Create GlnCheckBox
            app.GlnCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlnCheckBox.Enable = 'off';
            app.GlnCheckBox.Text = 'Gln';
            app.GlnCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GlnCheckBox.Position = [75 62 40 22];
            app.GlnCheckBox.Value = true;

            % Create GluCheckBox
            app.GluCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GluCheckBox.Enable = 'off';
            app.GluCheckBox.Text = 'Glu';
            app.GluCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GluCheckBox.Position = [134 62 40 22];
            app.GluCheckBox.Value = true;

            % Create GlyCheckBox
            app.GlyCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.GlyCheckBox.Enable = 'off';
            app.GlyCheckBox.Text = 'Gly';
            app.GlyCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.GlyCheckBox.Position = [194 62 40 22];

            % Create H2OCheckBox
            app.H2OCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.H2OCheckBox.Enable = 'off';
            app.H2OCheckBox.Text = 'H2O';
            app.H2OCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.H2OCheckBox.Position = [255 62 46 22];
            app.H2OCheckBox.Value = true;

            % Create InsCheckBox
            app.InsCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.InsCheckBox.Enable = 'off';
            app.InsCheckBox.Text = 'Ins';
            app.InsCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.InsCheckBox.Position = [316 62 38 22];
            app.InsCheckBox.Value = true;

            % Create LacCheckBox
            app.LacCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.LacCheckBox.Enable = 'off';
            app.LacCheckBox.Text = 'Lac';
            app.LacCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.LacCheckBox.Position = [373 62 42 22];
            app.LacCheckBox.Value = true;

            % Create NAACheckBox
            app.NAACheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAACheckBox.Enable = 'off';
            app.NAACheckBox.Text = 'NAA';
            app.NAACheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.NAACheckBox.Position = [427 62 46 22];
            app.NAACheckBox.Value = true;

            % Create NAAGCheckBox
            app.NAAGCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.NAAGCheckBox.Enable = 'off';
            app.NAAGCheckBox.Text = 'NAAG';
            app.NAAGCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.NAAGCheckBox.Position = [486 62 55 22];
            app.NAAGCheckBox.Value = true;

            % Create PChCheckBox
            app.PChCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PChCheckBox.Enable = 'off';
            app.PChCheckBox.Text = 'PCh';
            app.PChCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.PChCheckBox.Position = [545 62 45 22];
            app.PChCheckBox.Value = true;

            % Create PCrCheckBox
            app.PCrCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PCrCheckBox.Enable = 'off';
            app.PCrCheckBox.Text = 'PCr';
            app.PCrCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.PCrCheckBox.Position = [604 62 42 22];
            app.PCrCheckBox.Value = true;

            % Create PECheckBox
            app.PECheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PECheckBox.Enable = 'off';
            app.PECheckBox.Text = 'PE';
            app.PECheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.PECheckBox.Position = [655 62 37 22];
            app.PECheckBox.Value = true;

            % Create PhenylCheckBox
            app.PhenylCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.PhenylCheckBox.Enable = 'off';
            app.PhenylCheckBox.Text = 'Phenyl';
            app.PhenylCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.PhenylCheckBox.Position = [15 34 58 22];

            % Create ScylloCheckBox
            app.ScylloCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.ScylloCheckBox.Enable = 'off';
            app.ScylloCheckBox.Text = 'Scyllo';
            app.ScylloCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.ScylloCheckBox.Position = [75 34 54 22];
            app.ScylloCheckBox.Value = true;

            % Create SerCheckBox
            app.SerCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.SerCheckBox.Enable = 'off';
            app.SerCheckBox.Text = 'Ser';
            app.SerCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.SerCheckBox.Position = [134 34 40 22];

            % Create TauCheckBox
            app.TauCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TauCheckBox.Enable = 'off';
            app.TauCheckBox.Text = 'Tau';
            app.TauCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.TauCheckBox.Position = [194 34 41 22];
            app.TauCheckBox.Value = true;

            % Create TyrosCheckBox
            app.TyrosCheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.TyrosCheckBox.Enable = 'off';
            app.TyrosCheckBox.Text = 'Tyros';
            app.TyrosCheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.TyrosCheckBox.Position = [255 34 50 22];

            % Create MM09CheckBox
            app.MM09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM09CheckBox.Enable = 'off';
            app.MM09CheckBox.Text = 'MM09';
            app.MM09CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.MM09CheckBox.Position = [15 7 56 22];
            app.MM09CheckBox.Value = true;

            % Create MM12CheckBox
            app.MM12CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM12CheckBox.Enable = 'off';
            app.MM12CheckBox.Text = 'MM12';
            app.MM12CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.MM12CheckBox.Position = [75 7 56 22];
            app.MM12CheckBox.Value = true;

            % Create MM14CheckBox
            app.MM14CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM14CheckBox.Enable = 'off';
            app.MM14CheckBox.Text = 'MM14';
            app.MM14CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.MM14CheckBox.Position = [134 7 56 22];
            app.MM14CheckBox.Value = true;

            % Create MM17CheckBox
            app.MM17CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM17CheckBox.Enable = 'off';
            app.MM17CheckBox.Text = 'MM17';
            app.MM17CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.MM17CheckBox.Position = [194 7 56 22];
            app.MM17CheckBox.Value = true;

            % Create MM20CheckBox
            app.MM20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.MM20CheckBox.Enable = 'off';
            app.MM20CheckBox.Text = 'MM20';
            app.MM20CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.MM20CheckBox.Position = [255 7 56 22];
            app.MM20CheckBox.Value = true;

            % Create Lip09CheckBox
            app.Lip09CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip09CheckBox.Enable = 'off';
            app.Lip09CheckBox.Text = 'Lip09';
            app.Lip09CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.Lip09CheckBox.Position = [316 7 52 22];
            app.Lip09CheckBox.Value = true;

            % Create Lip13CheckBox
            app.Lip13CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip13CheckBox.Enable = 'off';
            app.Lip13CheckBox.Text = 'Lip13';
            app.Lip13CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.Lip13CheckBox.Position = [373 7 52 22];
            app.Lip13CheckBox.Value = true;

            % Create Lip20CheckBox
            app.Lip20CheckBox = uicheckbox(app.SelectedMetabolitesPanel);
            app.Lip20CheckBox.Enable = 'off';
            app.Lip20CheckBox.Text = 'Lip20';
            app.Lip20CheckBox.FontColor = [0.0392 0.2784 0.4392];
            app.Lip20CheckBox.Position = [427 7 52 22];
            app.Lip20CheckBox.Value = true;

            % Create SpecifyMRSandStructuralImagingFilesPanel
            app.SpecifyMRSandStructuralImagingFilesPanel = uipanel(app.InteractiveOspreyJobmCreatorUIFigure);
            app.SpecifyMRSandStructuralImagingFilesPanel.ForegroundColor = [0.0392 0.2784 0.4392];
            app.SpecifyMRSandStructuralImagingFilesPanel.Title = '3. Specify MRS and Structural Imaging Files';
            app.SpecifyMRSandStructuralImagingFilesPanel.BackgroundColor = [1 1 1];
            app.SpecifyMRSandStructuralImagingFilesPanel.FontWeight = 'bold';
            app.SpecifyMRSandStructuralImagingFilesPanel.FontSize = 15;
            app.SpecifyMRSandStructuralImagingFilesPanel.Position = [365 377 350 249];

            % Create MRSDataButton
            app.MRSDataButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.MRSDataButton.ButtonPushedFcn = createCallbackFcn(app, @MRSDataButtonPushed, true);
            app.MRSDataButton.BackgroundColor = [1 1 1];
            app.MRSDataButton.FontWeight = 'bold';
            app.MRSDataButton.FontColor = [0.0392 0.2784 0.4392];
            app.MRSDataButton.Position = [10 151 118 23];
            app.MRSDataButton.Text = 'MRS Data';

            % Create H2OReferenceButton
            app.H2OReferenceButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.H2OReferenceButton.ButtonPushedFcn = createCallbackFcn(app, @H2OReferenceButtonPushed, true);
            app.H2OReferenceButton.BackgroundColor = [1 1 1];
            app.H2OReferenceButton.FontWeight = 'bold';
            app.H2OReferenceButton.FontColor = [0.0392 0.2784 0.4392];
            app.H2OReferenceButton.Enable = 'off';
            app.H2OReferenceButton.Position = [10 119 118 23];
            app.H2OReferenceButton.Text = 'H2O Reference';

            % Create H2OShortTEButton
            app.H2OShortTEButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.H2OShortTEButton.ButtonPushedFcn = createCallbackFcn(app, @H2OShortTEButtonPushed, true);
            app.H2OShortTEButton.BackgroundColor = [1 1 1];
            app.H2OShortTEButton.FontWeight = 'bold';
            app.H2OShortTEButton.FontColor = [0.0392 0.2784 0.4392];
            app.H2OShortTEButton.Enable = 'off';
            app.H2OShortTEButton.Position = [10 84 118 23];
            app.H2OShortTEButton.Text = 'H2O Short TE';

            % Create MetaboliteNulledButton
            app.MetaboliteNulledButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.MetaboliteNulledButton.ButtonPushedFcn = createCallbackFcn(app, @MetaboliteNulledButtonPushed, true);
            app.MetaboliteNulledButton.BackgroundColor = [1 1 1];
            app.MetaboliteNulledButton.FontWeight = 'bold';
            app.MetaboliteNulledButton.FontColor = [0.0392 0.2784 0.4392];
            app.MetaboliteNulledButton.Enable = 'off';
            app.MetaboliteNulledButton.Position = [10 52 118 23];
            app.MetaboliteNulledButton.Text = 'Metabolite-Nulled';

            % Create T1DataniftiniiButton
            app.T1DataniftiniiButton = uibutton(app.SpecifyMRSandStructuralImagingFilesPanel, 'push');
            app.T1DataniftiniiButton.ButtonPushedFcn = createCallbackFcn(app, @T1DataniftiniiButtonPushed, true);
            app.T1DataniftiniiButton.BackgroundColor = [1 1 1];
            app.T1DataniftiniiButton.FontWeight = 'bold';
            app.T1DataniftiniiButton.FontColor = [0.0392 0.2784 0.4392];
            app.T1DataniftiniiButton.Enable = 'off';
            app.T1DataniftiniiButton.Position = [11 18 117 23];
            app.T1DataniftiniiButton.Text = 'T1 Data (nifti *.nii)';

            % Create NumberofdatasetsEditFieldLabel
            app.NumberofdatasetsEditFieldLabel = uilabel(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.NumberofdatasetsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofdatasetsEditFieldLabel.FontColor = [0.0392 0.2784 0.4392];
            app.NumberofdatasetsEditFieldLabel.Position = [18 191 112 22];
            app.NumberofdatasetsEditFieldLabel.Text = 'Number of datasets';

            % Create NumberofdatasetsEditField
            app.NumberofdatasetsEditField = uieditfield(app.SpecifyMRSandStructuralImagingFilesPanel, 'numeric');
            app.NumberofdatasetsEditField.Limits = [1 Inf];
            app.NumberofdatasetsEditField.ValueDisplayFormat = '%.0f';
            app.NumberofdatasetsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofdatasetsEditFieldValueChanged, true);
            app.NumberofdatasetsEditField.HorizontalAlignment = 'center';
            app.NumberofdatasetsEditField.FontColor = [0.0392 0.2784 0.4392];
            app.NumberofdatasetsEditField.Position = [145 191 100 22];
            app.NumberofdatasetsEditField.Value = 1;

            % Create MRSdataDropDown
            app.MRSdataDropDown = uidropdown(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.MRSdataDropDown.Items = {};
            app.MRSdataDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.MRSdataDropDown.BackgroundColor = [1 1 1];
            app.MRSdataDropDown.Position = [145 151 182 22];
            app.MRSdataDropDown.Value = {};

            % Create H2OReferenceDropDown
            app.H2OReferenceDropDown = uidropdown(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.H2OReferenceDropDown.Items = {};
            app.H2OReferenceDropDown.Enable = 'off';
            app.H2OReferenceDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.H2OReferenceDropDown.BackgroundColor = [1 1 1];
            app.H2OReferenceDropDown.Position = [145 119 182 22];
            app.H2OReferenceDropDown.Value = {};

            % Create H2OShortTEDropDown
            app.H2OShortTEDropDown = uidropdown(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.H2OShortTEDropDown.Items = {};
            app.H2OShortTEDropDown.Enable = 'off';
            app.H2OShortTEDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.H2OShortTEDropDown.BackgroundColor = [1 1 1];
            app.H2OShortTEDropDown.Position = [145 84 182 22];
            app.H2OShortTEDropDown.Value = {};

            % Create MetaboliteNulledDropDown
            app.MetaboliteNulledDropDown = uidropdown(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.MetaboliteNulledDropDown.Items = {};
            app.MetaboliteNulledDropDown.Enable = 'off';
            app.MetaboliteNulledDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.MetaboliteNulledDropDown.BackgroundColor = [1 1 1];
            app.MetaboliteNulledDropDown.Position = [145 52 182 22];
            app.MetaboliteNulledDropDown.Value = {};

            % Create T1DataniftiniiDropDown
            app.T1DataniftiniiDropDown = uidropdown(app.SpecifyMRSandStructuralImagingFilesPanel);
            app.T1DataniftiniiDropDown.Items = {};
            app.T1DataniftiniiDropDown.Enable = 'off';
            app.T1DataniftiniiDropDown.FontColor = [0.0392 0.2784 0.4392];
            app.T1DataniftiniiDropDown.BackgroundColor = [1 1 1];
            app.T1DataniftiniiDropDown.Position = [145 17 182 22];
            app.T1DataniftiniiDropDown.Value = {};

            % Create SpecifyOutputFolderPanel
            app.SpecifyOutputFolderPanel = uipanel(app.InteractiveOspreyJobmCreatorUIFigure);
            app.SpecifyOutputFolderPanel.ForegroundColor = [0.0392 0.2784 0.4392];
            app.SpecifyOutputFolderPanel.Title = '4. Specify Output Folder';
            app.SpecifyOutputFolderPanel.BackgroundColor = [1 1 1];
            app.SpecifyOutputFolderPanel.FontWeight = 'bold';
            app.SpecifyOutputFolderPanel.FontSize = 15;
            app.SpecifyOutputFolderPanel.Position = [365 257 350 107];

            % Create OutputFolderEditField
            app.OutputFolderEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.OutputFolderEditField.Editable = 'off';
            app.OutputFolderEditField.Position = [144 48 183 22];

            % Create OutputFolderButton
            app.OutputFolderButton = uibutton(app.SpecifyOutputFolderPanel, 'push');
            app.OutputFolderButton.ButtonPushedFcn = createCallbackFcn(app, @OutputFolderButtonPushed, true);
            app.OutputFolderButton.BackgroundColor = [1 1 1];
            app.OutputFolderButton.FontWeight = 'bold';
            app.OutputFolderButton.FontColor = [0.0392 0.2784 0.4392];
            app.OutputFolderButton.Position = [9 47 119 23];
            app.OutputFolderButton.Text = 'Output Folder';

            % Create JobNameEditFieldLabel
            app.JobNameEditFieldLabel = uilabel(app.SpecifyOutputFolderPanel);
            app.JobNameEditFieldLabel.HorizontalAlignment = 'right';
            app.JobNameEditFieldLabel.FontColor = [0.0392 0.2784 0.4392];
            app.JobNameEditFieldLabel.Position = [18 12 61 22];
            app.JobNameEditFieldLabel.Text = 'Job Name';

            % Create JobNameEditField
            app.JobNameEditField = uieditfield(app.SpecifyOutputFolderPanel, 'text');
            app.JobNameEditField.FontColor = [0.0392 0.2784 0.4392];
            app.JobNameEditField.Position = [144 12 183 22];

            % Create ExitButton
            app.ExitButton = uibutton(app.InteractiveOspreyJobmCreatorUIFigure, 'push');
            app.ExitButton.ButtonPushedFcn = createCallbackFcn(app, @ExitButtonPushed, true);
            app.ExitButton.BackgroundColor = [1 1 1];
            app.ExitButton.FontSize = 15;
            app.ExitButton.FontWeight = 'bold';
            app.ExitButton.FontColor = [0.0392 0.2784 0.4392];
            app.ExitButton.Position = [544 167 165 38];
            app.ExitButton.Text = 'Exit';

            % Create RunJobButton
            app.RunJobButton = uibutton(app.InteractiveOspreyJobmCreatorUIFigure, 'push');
            app.RunJobButton.ButtonPushedFcn = createCallbackFcn(app, @RunJobButtonPushed, true);
            app.RunJobButton.BackgroundColor = [1 1 1];
            app.RunJobButton.FontSize = 15;
            app.RunJobButton.FontWeight = 'bold';
            app.RunJobButton.FontColor = [0.0392 0.2784 0.4392];
            app.RunJobButton.Position = [365 214 168 38];
            app.RunJobButton.Text = 'Run Job';

            % Create LoadJobFileButton
            app.LoadJobFileButton = uibutton(app.InteractiveOspreyJobmCreatorUIFigure, 'push');
            app.LoadJobFileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadJobFileButtonPushed, true);
            app.LoadJobFileButton.BackgroundColor = [1 1 1];
            app.LoadJobFileButton.FontSize = 15;
            app.LoadJobFileButton.FontWeight = 'bold';
            app.LoadJobFileButton.FontColor = [0.0392 0.2784 0.4392];
            app.LoadJobFileButton.Position = [543 214 169 38];
            app.LoadJobFileButton.Text = 'Load JobFile';

            % Create LoadMRSContButton
            app.LoadMRSContButton = uibutton(app.InteractiveOspreyJobmCreatorUIFigure, 'push');
            app.LoadMRSContButton.ButtonPushedFcn = createCallbackFcn(app, @LoadMRSContButtonPushed, true);
            app.LoadMRSContButton.BackgroundColor = [1 1 1];
            app.LoadMRSContButton.FontSize = 15;
            app.LoadMRSContButton.FontWeight = 'bold';
            app.LoadMRSContButton.FontColor = [0.0392 0.2784 0.4392];
            app.LoadMRSContButton.Position = [366 167 167 38];
            app.LoadMRSContButton.Text = 'Load MRSCont';

            % Show the figure after all components are created
            app.InteractiveOspreyJobmCreatorUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CreateOspreyJob_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.InteractiveOspreyJobmCreatorUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.InteractiveOspreyJobmCreatorUIFigure)
        end
    end
end