function jobm = osp_create_job_file(app)
%%  osp_create_job_file
%   This function creates an Osprey jobFile based onn the entries made in the GUI.
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

outputFolder = app.OutputFolderEditField.Value;

if ~exist(outputFolder)
    mkdir(outputFolder)
end

if isempty(outputFolder)
    outputFolder = cd;
end

JobName = app.JobNameEditField.Value;

if isempty(JobName)
    JobName = 'jobGUIcreated';
end

mfile = fullfile(outputFolder,[JobName '.json']);

fid = fopen(mfile,'w');
fprintf(fid,'%s','{');
fprintf(fid,'\n\t%s',['"seqType": "' app.SequenceTypeDropDown.Value '",']);

if strcmp(app.SequenceTypeDropDown.Value,'MEGA')
    fprintf(fid,'\n\t%s',['"editTarget": ["' app.EditingTargetsDropDown.Value '"],']);
end
if strcmp(app.SequenceTypeDropDown.Value,'HERMES') || strcmp(app.SequenceTypeDropDown.Value,'HERCULES')
    if strcmp(app.EditingTargetsDropDown.Value, 'GABA, GSH')
        fprintf(fid,'\n\t%s',['"editTarget": ["GABA","GSH"],']);
    end
end

fprintf(fid,'\n\t%s',['"dataScenario": "' app.DataScenarioDropDown.Value '",']);
fprintf(fid,'\n\t%s',['"MM3coModel": "' app.MM3coDropDown.Value '",']);
fprintf(fid,'\n\t%s',['"FWHMMM3co": "' num2str(app.FWHMMM3coEditField.Value) '",']);
fprintf(fid,'\n\t%s',['"SpecReg": "' app.SpectralRegistrationDropDown.Value '",']);
fprintf(fid,'\n\t%s',['"SubSpecAlignment": "' app.SubspectraAligmentDropDown.Value '",']);
fprintf(fid,'\n\t%s',['"UnstableWater": "' num2str(app.unstablewaterCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"saveLCM": "' num2str(app.SaveLCMCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"savejMRUI": "' num2str(app.SaveJMRUICheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"saveVendor": "' num2str(app.SaveVendorCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"saveNII": "' num2str(app.SaveNIfTiCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"savePDF": "' num2str(app.SavePDFCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"method": "' app.FittingAlgorithmDropDown.Value '",']);
fprintf(fid,'\n\t%s',['"ECCmetab": "' num2str(app.ECCmetaboliteMRSCheckBox.Value) '",']);
fprintf(fid,'\n\t%s',['"ECCmm": "' num2str(app.ECCmmMRSCheckBox.Value) '",']);


switch app.IncludedMetabolitesDropDown.Value 
    case 'custom'
        first = 1;
        fprintf(fid,'\n\t%s','"includeMetabs": [');
        if app.AlaCheckBox.Value, fprintf(fid,'%s','"Ala"'), first = 0;end
        if app.AscCheckBox.Value, if first, fprintf(fid,'%s','"Asc"'), first = 0;else fprintf(fid,'%s',',"Asc"'); end; end
        if app.AspCheckBox.Value, if first, fprintf(fid,'%s','"Asp"'), first = 0;else fprintf(fid,'%s',',"Asp"'); end; end
        if app.bHBCheckBox.Value, if first, fprintf(fid,'%s','"bHB"'), first = 0;else fprintf(fid,'%s',',"bHB"'); end; end
        if app.bHGCheckBox.Value, if first, fprintf(fid,'%s','"bHG"'), first = 0;else fprintf(fid,'%s',',"bHG"'); end; end
        if app.CitCheckBox.Value, if first, fprintf(fid,'%s','"Cit"'), first = 0;else fprintf(fid,'%s',',"Cit"'); end; end
        if app.CrCheckBox.Value, if first, fprintf(fid,'%s','"Cr"'), first = 0;else fprintf(fid,'%s',',"Cr"'); end; end
        if app.CystatCheckBox.Value, fprintf(fid,'%s',',"Cystat"'); end
        if app.CrCH2CheckBox.Value, fprintf(fid,'%s',',"CrCH2"'); end
        if app.EtOHCheckBox.Value, fprintf(fid,'%s',',"EtOH"'); end
        if app.GABACheckBox.Value, fprintf(fid,'%s',',"GABA"'); end
        if app.GPCCheckBox.Value, fprintf(fid,'%s',',"GPC"'); end
        if app.GSHCheckBox.Value, fprintf(fid,'%s',',"GSH"'); end
        if app.GlcCheckBox.Value, fprintf(fid,'%s',',"Glc"'); end
        if app.GlnCheckBox.Value, fprintf(fid,'%s',',"Gln"'); end
        if app.GluCheckBox.Value, fprintf(fid,'%s',',"Glu"'); end
        if app.GlyCheckBox.Value, fprintf(fid,'%s',',"Gly"'); end
        if app.H2OCheckBox.Value, fprintf(fid,'%s',',"H2O"'); end
        if app.InsCheckBox.Value, fprintf(fid,'%s',',"mI"'); end
        if app.LacCheckBox.Value, fprintf(fid,'%s',',"Lac"'); end
        if app.NAACheckBox.Value, fprintf(fid,'%s',',"NAA"'); end
        if app.NAAGCheckBox.Value, fprintf(fid,'%s',',"NAAG"'); end
        if app.PChCheckBox.Value, fprintf(fid,'%s',',"PCh"'); end
        if app.PCrCheckBox.Value, fprintf(fid,'%s',',"PCr"'); end
        if app.PECheckBox.Value, fprintf(fid,'%s',',"PE"'); end
        if app.ScylloCheckBox.Value, fprintf(fid,'%s',',"sI"'); end
        if app.SerCheckBox.Value, fprintf(fid,'%s',',"Ser"'); end
        if app.TauCheckBox.Value, fprintf(fid,'%s',',"Tau"'); end
        if app.TyrosCheckBox.Value, fprintf(fid,'%s',',"Tyros",'); end
        if app.MM09CheckBox.Value, fprintf(fid,'\n\t\t\t\t\t\t\t%s',',"MM09"'); end
        if app.MM12CheckBox.Value, fprintf(fid,'%s',',"MM12"'); end
        if app.MM14CheckBox.Value, fprintf(fid,'%s',',"MM14"'); end
        if app.MM17CheckBox.Value, fprintf(fid,'%s',',"MM17"'); end
        if app.MM20CheckBox.Value, fprintf(fid,'%s',',"MM20"'); end
        if app.Lip09CheckBox.Value, fprintf(fid,'%s',',"Lip09"'); end
        if app.Lip13CheckBox.Value, fprintf(fid,'%s',',"Lip13"'); end
        if app.Lip20CheckBox.Value, fprintf(fid,'%s',',"Lip20"'); end
        if app.MM37CheckBox.Value, fprintf(fid,'%s',',"MM37"'); end
        if app.MM38CheckBox.Value, fprintf(fid,'%s',',"MM38"'); end
        if app.MM40CheckBox.Value, fprintf(fid,'%s',',"MM40"'); end
        if app.MM42CheckBox.Value, fprintf(fid,'%s',',"MM42"'); end
        if app.MMexpCheckBox.Value, fprintf(fid,'%s',',"MMexp"'); end
        if app.MM_PCCCheckBox.Value, fprintf(fid,'%s',',"MM_PRESS_PCC"'); end
        if app.MM_CSOCheckBox.Value, fprintf(fid,'%s',',"MM_PRESS_CSO"'); end
        fprintf(fid,'%s','],');
    
    otherwise
        fprintf(fid,'\n\t%s',['"includeMetabs": ["' app.IncludedMetabolitesDropDown.Value '"],']);
end

fprintf(fid,'\n\t%s',['"style": "' app.FittingStyleDropDown.Value '",']);
range = strsplit(app.MRSFitRangeppmEditField.Value);
fprintf(fid,'\n\t%s',['"lolim_range": "' num2str(range{1}) '",']);
fprintf(fid,'\n\t%s',['"uplim_range": "' num2str(range{2}) '",']);
rangew = strsplit(app.WaterFitRangeppmEditField.Value);
fprintf(fid,'\n\t%s',['"lolim_rangew": "' num2str(num2str(rangew{1})) '",']);
fprintf(fid,'\n\t%s',['"uplim_rangew": "' num2str(num2str(rangew{2})) '",']);
fprintf(fid,'\n\t%s',['"bLineKnotSpace": "' num2str(app.BaselineknotspacingppmEditField.Value) '",']);
fprintf(fid,'\n\t%s',['"fitMM": "' num2str(app.AddMMandLipbasisfunctionstofitCheckBox.Value) '",']);

if ~isempty(app.BasisSetEditField.Value)
    fprintf(fid,'\n\t%s',['"basisSet": "' app.BasisSetEditField.Value '",']);
end

if isempty(app.MRSDataText.Value{1})
    error('A MRS data file should be specified')
else
    fprintf(fid,'\n\t%s',['"files": ["' app.MRSDataText.Value{1} '"']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',',');
        fprintf(fid,'\n\t%s',['"' app.MRSDataText.Value{i} '"']);
    end
    fprintf(fid,'%s','],');
end

if ~isempty(app.H2OReferenceText.Value{1})
    fprintf(fid,'\n\t%s',['"files_ref": ["' app.H2OReferenceText.Value{1} '"']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',',');
        fprintf(fid,'\n\t%s',['"' app.H2OReferenceText.Value{i} '"']);
    end
    fprintf(fid,'%s','],');
end

if ~isempty(app.H2OShortTEText.Value{1})
    fprintf(fid,'\n\t%s',['"files_w": ["' app.H2OShortTEText.Value{1} '"']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',',');
        fprintf(fid,'\n\t%s',['"' app.H2OShortTEText.Value{i} '"']);
    end
    fprintf(fid,'%s','],');
end

if ~isempty(app.MetaboliteNulledText.Value{1})
    fprintf(fid,'\n\t%s',['"files_mm": ["' app.MetaboliteNulledText.Value{1} '"']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',',');
        fprintf(fid,'\n\t%s',['"' app.MetaboliteNulledText.Value{i} '"']);
    end
    fprintf(fid,'%s','],');
end

if ~isempty(app.T1DataText.Value{1})
    fprintf(fid,'\n\t%s',['"files_nii": ["' app.T1DataText.Value{1} '"']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',',');
        fprintf(fid,'\n\t%s',['"' app.T1DataText.Value{i} '"']);
    end
    fprintf(fid,'%s','],');
end

fprintf(fid,'\n\t%s',['"file_stat": ["' app.StatcsvEditField.Value '"],']);
fprintf(fid,'\n\t%s',['"outputFolder": ["' app.OutputFolderEditField.Value '"]']);
fprintf(fid,'\n%s','}');
fclose(fid);

jobm = mfile;