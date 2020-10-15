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

if isempty(outputFolder)
    outputFolder = cd;
end

JobName = app.JobNameEditField.Value;

if isempty(JobName)
    JobName = 'jobGUIcreated';
end

mfile = fullfile(outputFolder,[JobName '.m']);

fid = fopen(mfile,'w');
fprintf(fid,'%s','%%% 0. CREDITS %%%');
fprintf(fid,'\n%s','% This is a GUI generated Osprey jobFile. The code was kindly shared by Dr. Peter Van Schuerbeek (UZ Brussel)');
fprintf(fid,'\n%s','');
fprintf(fid,'\n%s','%%% 1. SPECIFY SEQUENCE INFORMATION %%%');
fprintf(fid,'\n%s',['seqType = ''' app.SequenceTypeDropDown.Value ''';']);
fprintf(fid,'\n%s',['editTarget = {''' app.EditingTargetsDropDown.Value '''};']);
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 2. SPECIFY DATA HANDLING AND MODELING OPTIONS %%%');
fprintf(fid,'\n%s',['opts.saveLCM = ' num2str(app.SaveLCMCheckBox.Value) ';']);
fprintf(fid,'\n%s',['opts.savejMRUI = ' num2str(app.SaveJMRUICheckBox.Value) ';']);
fprintf(fid,'\n%s',['opts.saveVendor = ' num2str(app.SaveVendorCheckBox.Value) ';']);
fprintf(fid,'\n%s',['opts.protocol = ''' app.MRSProtocolDropDown.Value ''';']);

switch lower(app.SpectralregistrationDropDown.Value)
    case 'robust'
        fprintf(fid,'\n%s',['opts.SpecReg = ''RobSpecReg'';']);
    case 'restricted'
        fprintf(fid,'\n%s',['opts.SpecReg = ''RestrSpecReg'';']);
    case 'none'
        fprintf(fid,'\n%s',['opts.SpecReg = ''none'';']);
end

fprintf(fid,'\n%s',['opts.denoising = ''' app.DenoisingDropDown.Value ''';']);

fprintf(fid,'\n%s',['opts.fit.method = ''' app.FittingAlgorithmDropDown.Value ''';']);

switch app.IncludedMetabolitesDropDown.Value 
    case 'custom'
        fprintf(fid,'\n%s','opts.fit.includeMetabs = {');
        if app.AlaCheckBox.Value, fprintf(fid,'%s',' ''Ala'''); end
        if app.AscCheckBox.Value, fprintf(fid,'%s',' ''Asc'''); end
        if app.AspCheckBox.Value, fprintf(fid,'%s',' ''Asp'''); end
        if app.bHBCheckBox.Value, fprintf(fid,'%s',' ''bHB'''); end
        if app.bHGCheckBox.Value, fprintf(fid,'%s',' ''bHG'''); end
        if app.CitCheckBox.Value, fprintf(fid,'%s',' ''Cit'''); end
        if app.CrCheckBox.Value, fprintf(fid,'%s',' ''Cr'''); end
        if app.CrCH2CheckBox.Value, fprintf(fid,'%s',' ''CrCH2'''); end
        if app.EtOHCheckBox.Value, fprintf(fid,'%s',' ''EtOH'''); end
        if app.GABACheckBox.Value, fprintf(fid,'%s',' ''GABA'''); end
        if app.GPCCheckBox.Value, fprintf(fid,'%s',' ''GPC'''); end
        if app.GSHCheckBox.Value, fprintf(fid,'%s',' ''GSH'''); end
        if app.GlcCheckBox.Value, fprintf(fid,'%s',' ''Glc'''); end
        if app.GlnCheckBox.Value, fprintf(fid,'%s',' ''Gln'''); end
        if app.GluCheckBox.Value, fprintf(fid,'%s',' ''Glu'''); end
        if app.GlyCheckBox.Value, fprintf(fid,'%s',' ''Gly'''); end
        if app.H2OCheckBox.Value, fprintf(fid,'%s',' ''H2O'''); end
        if app.InsCheckBox.Value, fprintf(fid,'%s',' ''Ins'''); end
        if app.LacCheckBox.Value, fprintf(fid,'%s',' ''Lac'''); end
        if app.NAACheckBox.Value, fprintf(fid,'%s',' ''NAA'''); end
        if app.NAAGCheckBox.Value, fprintf(fid,'%s',' ''NAAG'''); end
        if app.PChCheckBox.Value, fprintf(fid,'%s',' ''PCh'''); end
        if app.PCrCheckBox.Value, fprintf(fid,'%s',' ''PCr'''); end
        if app.PECheckBox.Value, fprintf(fid,'%s',' ''PE'''); end
        if app.PhenylCheckBox.Value, fprintf(fid,'%s',' ''Phenyl'''); end
        if app.ScylloCheckBox.Value, fprintf(fid,'%s',' ''Scyllo'''); end
        if app.SerCheckBox.Value, fprintf(fid,'%s',' ''Ser'''); end
        if app.TauCheckBox.Value, fprintf(fid,'%s',' ''Tau'''); end
        if app.TyrosCheckBox.Value, fprintf(fid,'%s',' ''Tyros'''); end
        if app.MM09CheckBox.Value, fprintf(fid,'%s',' ''MM09'''); end
        if app.MM12CheckBox.Value, fprintf(fid,'%s',' ''MM12'''); end
        if app.MM14CheckBox.Value, fprintf(fid,'%s',' ''MM14'''); end
        if app.MM17CheckBox.Value, fprintf(fid,'%s',' ''MM17'''); end
        if app.MM20CheckBox.Value, fprintf(fid,'%s',' ''MM20'''); end
        if app.Lip09CheckBox.Value, fprintf(fid,'%s',' ''Lip09'''); end
        if app.Lip13CheckBox.Value, fprintf(fid,'%s',' ''Lip13'''); end
        if app.Lip20CheckBox.Value, fprintf(fid,'%s',' ''Lip20'''); end
        fprintf(fid,'%s','};');
    
    otherwise
        fprintf(fid,'\n%s',['opts.fit.includeMetabs = {''' app.IncludedMetabolitesDropDown.Value '''};']);
end

fprintf(fid,'\n%s',['opts.fit.style = ''' app.FittingStyleDropDown.Value ''';']);
fprintf(fid,'\n%s',['opts.fit.reffreq = ' num2str(app.ReferencefrequencyppmEditField.Value) ';']);
fprintf(fid,'\n%s',['opts.fit.range = [' app.MRSFitRangeppmEditField.Value '];']);
fprintf(fid,'\n%s',['opts.fit.rangeWater = [' app.WaterFitRangeppmEditField.Value '];']);
fprintf(fid,'\n%s',['opts.fit.bLineKnotSpace = ' num2str(app.BaselineknotspacingppmEditField.Value) ';']);
fprintf(fid,'\n%s',['opts.fit.fitMM = ' num2str(app.AddMMandLipbasisfunctionstofitCheckBox.Value) ';']);
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 3. SPECIFY MRS DATA AND STRUCTURAL IMAGING FILES %%%');
if isempty(app.MRSdataDropDown.Items)
    error('A MRS data file should be specified')
else
    fprintf(fid,'\n%s',['files = {']);
    fprintf(fid,'%s',['''' app.MRSdataDropDown.Items{1} '''']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',[',...']);
        fprintf(fid,'\n\t\t%s',[' ' '''' app.MRSdataDropDown.Items{i} '''']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.H2OReferenceDropDown.Items)
    fprintf(fid,'\n%s',['files_ref = {};']);
else
    fprintf(fid,'\n%s',['files_ref = {']);
    fprintf(fid,'%s',['''' app.H2OReferenceDropDown.Items{1} '''']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',[',...']);
        fprintf(fid,'\n\t\t%s',[' ' '''' app.H2OReferenceDropDown.Items{i} '''']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.H2OShortTEDropDown.Items)
    fprintf(fid,'\n%s',['files_w = {};']);
else
    fprintf(fid,'\n%s',['files_w = {']);
    fprintf(fid,'%s',['''' app.H2OShortTEDropDown.Items{1} '''']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',[',...']);
        fprintf(fid,'\n\t\t%s',[' ' '''' app.H2OShortTEDropDown.Items{i} '''']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.MetaboliteNulledDropDown.Items)
    fprintf(fid,'\n%s',['files_mm = {};']);
else
    fprintf(fid,'\n%s',['files_mm = {']);
    fprintf(fid,'%s',['''' app.MetaboliteNulledDropDown.Items{1} '''']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',[',...']);
        fprintf(fid,'\n\t\t%s',[' ' '''' app.MetaboliteDropDown.Items{i} '''']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.T1DataniftiniiDropDown.Items)
    fprintf(fid,'\n%s',['files_nii = {};']);
else
    fprintf(fid,'\n%s',['files_nii = {']);
    fprintf(fid,'%s',['''' app.T1DataniftiniiDropDown.Items{1} '''']);
    for i=2:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',[',...']);
        fprintf(fid,'\n\t\t%s',[' ' '''' app.T1DataniftiniiDropDown.Items{i} '''']);
    end
    fprintf(fid,'%s',['};']);
end
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 4. SPECIFY OUTPUT FOLDER %%%');
fprintf(fid,'\n%s',['outputFolder = ''' app.OutputFolderEditField.Value ''';']);

fclose(fid);

jobm = mfile;