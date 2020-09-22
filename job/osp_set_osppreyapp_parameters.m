function osp_set_osppreyapp_parameters(mfile,mpath,app)

run(fullfile(mpath,mfile));

%%% 1. SPECIFY SEQUENCE INFORMATION %%%
app.SequenceTypeDropDown.Value = seqType;

switch seqType
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
            
app.EditingTargetsDropDown.Value = editTarget;

%%% 2. SPECIFY DATA HANDLING AND MODELING OPTIONS %%%
app.SaveLCMCheckBox.Value = opts.saveLCM;
app.SaveJMRUICheckBox.Value = opts.savejMRUI;
app.SaveVendorCheckBox.Value = opts.saveVendor;
if isfield(opts,'protocol')
    app.MRSProtocolDropDown.Value = opts.protocol;
else
    app.MRSProtocolDropDown.Value = 'Brain';
end

if isfield(opts,'SpecReg')
    switch opts.SpecReg
        case 'RobSpecReg'
            app.SpectralregistrationDropDown.Value = 'Robust';
        case 'RestrSpecReg'
            app.SpectralregistrationDropDown.Value = 'Restricted';
        case 'none' 
            app.SpectralregistrationDropDown.Value = 'Nonne';
    end
else
    if isfield(opts,'protocol')
        switch lower(opts.protocol)
            case 'braino pphantom'
                app.SpectralregistrationDropDown.Value = 'Restricted';
            otherwise
                app.SpectralregistrationDropDown.Value = 'Robust';
        end
    else
        app.SpectralregistrationDropDown.Value = 'Robust';
    end
end

app.FittingAlgorithmDropDown.Value = opts.fit.method;

switch opts.fit.includeMetabs{1}
    case 'custom'
        app.IncludedMetabolitesDropDown.Value = 'custom';
            
        if any(strcmp(opts.fit.includeMetabs,'Ala')); app.AlaCheckBox.Value = true; else app.AlaCheckBox.Value = false; end
        app.AlaCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Asc')); app.AscCheckBox.Value = true; else app.AscCheckBox.Value = false; end
        app.AscCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Asp')); app.AspCheckBox.Value = true; else app.AspCheckBox.Value = false; end
        app.AspCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'bHB')); app.bHBCheckBox.Value = true; else app.bHBCheckBox.Value = false; end
        app.bHBCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'bHG')); app.bHGCheckBox.Value = true; else app.bHGCheckBox.Value = false; end
        app.bHGCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Cit')); app.CitCheckBox.Value = true; else app.CitCheckBox.Value = false; end
        app.CitCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Cr')); app.CrCheckBox.Value = true; else app.CrCheckBox.Value = false; end
        app.CrCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'CrCH2')); app.CrCH2CheckBox.Value = true; else app.CrCH2CheckBox.Value = false; end
        app.CrCH2CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'EtOH')); app.EtOHCheckBox.Value = true; else app.EtOHCheckBox.Value = false; end
        app.EtOHCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'GABA')); app.GABACheckBox.Value = true; else app.GABACheckBox.Value = false; end
        app.GABACheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'GPC')); app.GPCCheckBox.Value = true; else app.GPCCheckBox.Value = false; end
        app.GPCCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'GSH')); app.GSHCheckBox.Value = true; else app.GSHCheckBox.Value = false; end
        app.GSHCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Glc')); app.GlcCheckBox.Value = true; else app.GlcCheckBox.Value = false; end
        app.GlcCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Gln')); app.GlnCheckBox.Value = true; else app.GlnCheckBox.Value = false; end
        app.GlnCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Glu')); app.GluCheckBox.Value = true; else app.GluCheckBox.Value = false; end
        app.GluCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Gly')); app.GlyCheckBox.Value = true; else app.GlyCheckBox.Value = false; end
        app.GlyCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'H2O')); app.H2OCheckBox.Value = true; else app.H2OCheckBox.Value = false; end
        app.H2OCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Ins')); app.InsCheckBox.Value = true; else app.InsCheckBox.Value = false; end
        app.InsCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Lac')); app.LacCheckBox.Value = true; else app.LacCheckBox.Value = false; end
        app.LacCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'NAA')); app.NAACheckBox.Value = true; else app.NAACheckBox.Value = false; end
        app.NAACheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'NAAG')); app.NAAGCheckBox.Value = true; else app.NAAGCheckBox.Value = false; end
        app.NAAGCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'PCh')); app.PChCheckBox.Value = true; else app.PChCheckBox.Value = false; end
        app.PChCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'PCr')); app.PCrCheckBox.Value = true; else app.PCrCheckBox.Value = false; end
        app.PCrCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'PE')); app.PECheckBox.Value = true; else app.PECheckBox.Value = false; end
        app.PECheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Phenyl')); app.PhenylCheckBox.Value = true; else app.PhenylCheckBox.Value = false; end
        app.PhenylCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Scyllo')); app.ScylloCheckBox.Value = true; else app.ScylloCheckBox.Value = false; end
        app.ScylloCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Ser')); app.SerCheckBox.Value = true; else app.SerCheckBox.Value = false; end
        app.SerCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Tau')); app.TauCheckBox.Value = true; else app.TauCheckBox.Value = false; end
        app.TauCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Tyros')); app.TyrosCheckBox.Value = true; else app.TyrosCheckBox.Value = false; end
        app.TyrosCheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'MM09')); app.MM09CheckBox.Value = true; else app.MM09CheckBox.Value = false; end
        app.MM09CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'MM12')); app.MM12CheckBox.Value = true; else app.MM12CheckBox.Value = false; end
        app.MM12CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'MM14')); app.MM14CheckBox.Value = true; else app.MM14CheckBox.Value = false; end
        app.MM14CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'MM17')); app.MM17CheckBox.Value = true; else app.MM17CheckBox.Value = false; end
        app.MM17CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'MM20')); app.MM20CheckBox.Value = true; else app.MM20CheckBox.Value = false; end
        app.MM20CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Lip09')); app.Lip09CheckBox.Value = true; else app.Lip09CheckBox.Value = false; end
        app.Lip09CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Lip13')); app.Lip13CheckBox.Value = true; else app.Lip13CheckBox.Value = false; end
        app.Lip13CheckBox.Enable = 'On';
        if any(strcmp(opts.fit.includeMetabs,'Lip20')); app.Lip20CheckBox.Value = true; else app.Lip20CheckBox.Value = false; end
        app.Lip20CheckBox.Enable = 'On';
    otherwise
        app.IncludedMetabolitesDropDown.Value = opts.fit.includeMetabs{1};
        
        osp_set_metabolites(app)
end

app.FittingStyleDropDown.Value = opts.fit.style;

if isfield(opts.fit,'reffreq')
    app.ReferencefrequencyppmEditField.Value = opts.fit.reffreq;
else
    if isfield(opts,'protocol')
        switch lower(opts.pprotocol)
            case 'braino phantom'
                app.ReferencefrequencyppmEditField.Value = 4.83;
            otherwise
                app.ReferencefrequencyppmEditField.Value = 4.68;
        end
    else
        app.ReferencefrequencyppmEditField.Value = 4.68;
    end
end

app.MRSFitRangeppmEditField.Value = [num2str(opts.fit.range(1),'%2.1f') ' ' num2str(opts.fit.range(2),'%2.1f')];
app.WaterFitRangeppmEditField.Value = [num2str(opts.fit.rangeWater(1),'%2.1f') ' ' num2str(opts.fit.rangeWater(2),'%2.1f')];
app.BaselineknotspacingppmEditField.Value = opts.fit.bLineKnotSpace;
app.AddMMandLipbasisfunctionstofitCheckBox.Value = opts.fit.fitMM;

%%% 3. SPECIFY MRS DATA AND STRUCTURAL IMAGING FILES %%%

app.NumberofdatasetsEditField.Value = numel(files);

if isempty(files)
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
else    
    app.MRSdataDropDown.Items = files;
    app.MRSdataDropDown.Value = files{1};
    
    [~,~,file_exten]=fileparts(files{1});
    
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
    
    app.H2OShortTEDropDown.Items = {};
    app.H2OShortTEDropDown.Value = {};
    app.H2OShortTEDropDown.Enable = 'On';
    app.H2OShortTEButton.Enable = 'On';
            
    app.MetaboliteNulledDropDown.Items = {};
    app.MetaboliteNulledDropDown.Value = {};
    app.MetaboliteNulledDropDown.Enable = 'On';
    app.MetaboliteNulledButton.Enable = 'On';
            
    app.T1DataniftiniiDropDown.Items = {};
    app.T1DataniftiniiDropDown.Value = {};
    app.T1DataniftiniiDropDown.Enable = 'On';
    app.T1DataniftiniiButton.Enable = 'On';
end
if isempty(files_ref)
    app.H2OReferenceDropDown.Items = {};
    app.H2OReferenceDropDown.Value = {};
else    
    app.H2OReferenceDropDown.Items = files;
    app.H2OReferenceDropDown.Value = files{1};
end
if isempty(files_w)
    app.H2OShortTEDropDown.Items = {};
    app.H2OShortTEDropDown.Value = {};
else    
    app.H2OShortTEDropDown.Items = files;
    app.H2OShortTEDropDown.Value = files{1};
end
if isempty(files_mm)
    app.MetaboliteNulledDropDown.Items = {};
    app.MetaboliteNulledDropDown.Value = {};
else    
    app.MetaboliteNulledDropDown.Items = files;
    app.MetaboliteNulledDropDown.Value = files{1};
end
if isempty(files_nii)
    app.T1DataniftiniiDropDown.Items = {};
    app.T1DataniftiniiDropDown.Value = {};
else    
    app.T1DataniftiniiDropDown.Items = files;
    app.T1DataniftiniiDropDown.Value = files{1};
end

%%% 4. SPECIFY OUTPUT FOLDER %%%
app.OutputFolderEditField.Value = outputFolder;

[~,file_name,file_exten]=fileparts(fullfile(mpath,mfile));

app.JobNameEditField.Value = file_name;

end