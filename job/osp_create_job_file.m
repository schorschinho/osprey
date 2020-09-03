function jobm = osp_create_job_file(app)

outputFolder = app.OutputFolderEditField.Value;

if isempty(outputFolder)
    outputFolder = cd;
end

JobName = app.JobNameEditField.Value;

mfile = fullfile(outputFolder,['ospreyJob_' JobName '.m']);

fid = fopen(mfile,'w');

fprintf(fid,'%s','%%% 1. SPECIFY SEQUENCE INFORMATION %%%');
fprintf(fid,'\n%s',['seqType = ''' app.SequenceTypeDropDown.Value ''';']);
fprintf(fid,'\n%s',['editTarget = {''' app.EditingTargetsDropDown.Value '''};']);
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 2. SPECIFY DATA HANDLING AND MODELING OPTIONS %%%');
fprintf(fid,'\n%s',['opts.saveLCM = ' num2str(app.SaveLCMCheckBox.Value) ';']);
fprintf(fid,'\n%s',['opts.savejMRUI = ' num2str(app.SaveJMRUICheckBox.Value) ';']);
fprintf(fid,'\n%s',['opts.saveVendor = ' num2str(app.SaveVendorCheckBox.Value) ';']);
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
fprintf(fid,'\n%s',['opts.fit.range = [' app.MRSFitRangeppmEditField.Value '];']);
fprintf(fid,'\n%s',['opts.fit.rangeWater = [' app.WaterFitRangeppmEditField.Value '];']);
fprintf(fid,'\n%s',['opts.fit.bLineKnotSpace = ' num2str(app.BaselineknotspacingppmEditField.Value) ';']);
fprintf(fid,'\n%s',['opts.fit.fitMM = ' num2str(app.AddMMandLipbasisfunctionstofitCheckBox.Value) ';']);
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 3. SPECIFY MRS DATA AND STRUCTURAL IMAGING FILES %%%');
if isempty(app.MRSDataText.Value{1})
    error('A MRS data file should be specified')
else
    fprintf(fid,'\n%s',['files = {']);
    for i=1:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',['''' app.MRSDataText.Value{i} ''',']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.H2OReferenceText.Value{1})
    fprintf(fid,'\n%s',['files_ref = {};']);
else
    fprintf(fid,'\n%s',['files_ref = {']);
    for i=1:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',['''' app.MRSDataText.Value{i} ''',']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.H2OShortTEText.Value{1})
    fprintf(fid,'\n%s',['files_w = {};']);
else
    fprintf(fid,'\n%s',['files_w = {']);
    for i=1:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',['''' app.H2OShortTEText.Value{i} ''',']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.MetaboliteNulledText.Value{1})
    fprintf(fid,'\n%s',['files_mm = {};']);
else
    fprintf(fid,'\n%s',['files_mm = {']);
    for i=1:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',['''' app.MetaboliteNulledText.Value{i} ''',']);
    end
    fprintf(fid,'%s',['};']);
end
if isempty(app.T1DataText.Value{1})
    fprintf(fid,'\n%s',['files_nii = {};']);
else
    fprintf(fid,'\n%s',['files_nii = {']);
    for i=1:app.NumberofdatasetsEditField.Value
        fprintf(fid,'%s',['''' app.T1DataText.Value{i} ''',']);
    end
    fprintf(fid,'%s',['};']);
end
fprintf(fid,'\n%s','');

fprintf(fid,'\n%s','%%% 4. SPECIFY OUTPUT FOLDER %%%');
fprintf(fid,'\n%s',['outputFolder = ''' app.OutputFolderEditField.Value ''';']);

fclose(fid);

jobm = mfile;