function MRSCont = osp_writelcm_userdefined_control(MRSCont, kk, which, LCMparam)
%% [MRSCont] = osp_writelcm_control(MRSCont, kk, which, LCMparam)
%   This function creates a LCModel compatible .control file which can be used
%   for LCModel batch processing. 
%
%   USAGE:
%       RF = osp_writelcm_control(MRSCont, kk, rr, which, LCMparam);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = Index for the kk-th dataset (optional. Default = 1)
%       which       = String for the spectrum to plot (optional)
%                   OPTIONS:    'A' (default)
%                               'B'
%                               'C'
%                               'D'
%                               'diff1'
%                               'diff2'
%                               'sum'
%      LCMparam     = struct with the control file parameters this struct created
%                     by the LCMcontrol.m script included in the Osprey process folder
%       
%                       
%   
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%      
%   ADHOC ADDITIONS    
%   Dr. Muhammad Saleh (University of Maryland Baltimore, 2022-05-04)
%   Allows user-defined destination folders and control parameters. You can run this function
%   independently after completing Osprey modules. 
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)


% Fall back to defaults if not provided
if nargin < 4
    error('ERROR: no input LCModel control struct specified.  Aborting!!');
    if nargin < 3
        which = 'A';
        if nargin < 2
            kk = 1;
            if nargin < 1
                error('ERROR: no input Osprey container specified.  Aborting!!');
            end
        end
    end
end
RF = LCMparam;

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'LCModelControlFiles');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

switch which
    case {'A','B','C','D'}
        subspec = ['outfile' which];
    case {'diff1'}
        subspec = 'outfileDiff1';
    case {'diff2'}
        subspec = 'outfileDiff2';
    case {'sum'}
        subspec = 'outfileSum';
end
name_raw = MRSCont.opts.fit.lcmodel.(subspec){kk};

% Retrieve the filenames of the LCModel .RAW files that were created at the
% end of OspreyProcess using osp_saveLCM.
if isfield(MRSCont.opts.fit.lcmodel, 'outfileRef')
    name_ref    = MRSCont.opts.fit.lcmodel.outfileRef{kk};
end
if isfield(MRSCont.opts.fit.lcmodel, 'outfileW')
    name_w      = MRSCont.opts.fit.lcmodel.outfileW{kk};
end


% If only water-suppressed data has been provided, set the LCModel control
% parameter for water-scaling to 'false'.
if ~exist('name_ref','var') && ~exist('name_w','var')
    LCMparam = osp_editControlParameters(LCMparam, 'dows', 'F');
else
    % If same-TE water reference data are provided, pass it on to LCModel, but
    % if ALSO a short-TE water reference has been provided, pass THAT ONE on.
    if exist('name_ref','var') && ~exist('name_w','var')
        selectedWaterRef = name_ref;
    elseif ~exist('name_ref','var') && exist('name_w','var')
        selectedWaterRef = name_w;
    elseif exist('name_ref','var') && exist('name_w','var')
        selectedWaterRef = name_w;
    end
end

% Add the necessary parameters
LCMparam = osp_editControlParameters(LCMparam, 'hzpppm', sprintf('%2.6e', MRSCont.processed.(which){kk}.txfrq/1e6));
LCMparam = osp_editControlParameters(LCMparam, 'deltat', sprintf('%2.6e', MRSCont.processed.(which){kk}.dwelltime));
LCMparam = osp_editControlParameters(LCMparam, 'nunfil', sprintf('%i', MRSCont.processed.(which){kk}.sz(1)));
if isfield(LCMparam, 'chcomb')
    LCMparam = osp_editControlParameters(LCMparam, 'ncombi', sprintf('%i', length(LCMparam.chcomb)));
end
if isfield(LCMparam, 'chomit')
    LCMparam = osp_editControlParameters(LCMparam, 'nomit', sprintf('%i', length(LCMparam.chomit)));
end

% Create a control file with the same name as the water-suppressed file.
[~, rawFilename, ~] = fileparts(name_raw);
%Modified to save within the subject folder -- 15Apr2022 mgs
controlFile = fullfile(strrep(name_raw, '.RAW', '.control'))
%Created a new other info -- 15Apr2022 MGSaleh
lcm_quant_File = fullfile(strrep(name_raw, '.RAW', ''));
% Create the LCModel output filenames (coord, print, ps, table, csv)
% Since the names are identical (and only the extension will vary), create
% the string without extension.
%Modified to save within the PSPath folder -- 15Apr2022 mgs
outputFilesCommon = fullfile(MRSCont.opts.pspath(kk,:), [MRSCont.opts.file_name(kk,:) '_' which]);

% Add individual LCM control parameters:
LCMparam = osp_editControlParameters(LCMparam, 'title', ['''' rawFilename '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filraw', ['''' name_raw '''']);     % path to the .RAW file produced previously
LCMparam = osp_editControlParameters(LCMparam, 'filh2o', ['''' name_w '''']);     % path to the .RAW file produced previously
LCMparam = osp_editControlParameters(LCMparam, 'filps',  ['''' strcat(outputFilesCommon, '.ps') '''']);
%Modified to save within the subject LCM folder -- 15Apr2022 mgs
LCMparam = osp_editControlParameters(LCMparam, 'filtab', ['''' strcat(lcm_quant_File, '.table') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filcsv', ['''' strcat(lcm_quant_File, '.csv') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filcoo', ['''' strcat(lcm_quant_File, '.coord') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filpri', ['''' strcat(lcm_quant_File, '.print') '''']);

if exist('name_ref','var')
    LCMparam = osp_editControlParameters(LCMparam, 'filh2o',  ['''' selectedWaterRef '''']);
end    

% Write every line into the control file
fid = fopen(controlFile,'w+');
fprintf(fid,' $LCMODL');
fprintf(fid,' \n OWNER=''Diagnostic Radiology Department, University of Maryland Baltimore'''); %Added the institute -- 15Apr2022 MGSaleh
% Create list of fieldnames and loop over them
allFields = fieldnames(LCMparam);
for ff = 1:length(allFields)
    % If the value is a cell, we need to loop over the elements
    if iscell(LCMparam.(allFields{ff}))
        for cc = 1:length(LCMparam.(allFields{ff}))
            fprintf(fid, '\n %s = %s', [allFields{ff} '(' num2str(cc) ')'], LCMparam.(allFields{ff}){cc});
        end
    else
        % If the value is simply a string, we can just write it.
        fprintf(fid, '\n %s = %s', allFields{ff}, LCMparam.(allFields{ff}));
    end
end
%General printouts -- 15Apr2022 MGSaleh
% fprintf(fid,'\n doecc=F');
fprintf(fid,'\n LTABLE=7'); %make table file
fprintf(fid,'\n ipage2 = 1');
fprintf(fid,'\n lprint = 6');
fprintf(fid,'\n Lcoord = 9');
fprintf(fid,'\n lcsv = 11');
fprintf(fid,'\n lps = 8');
fprintf(fid,'\n Neach = 50');
fprintf(fid,'\n ECHOT=%d',MRSCont.raw{1,kk}.te); 
switch which
    case {'diff1'}
        if MRSCont.flags.isHERMES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERMES_GABA.basis');
        elseif MRSCont.flags.isHERCULES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERC_GABA.basis');
        end
% %         fprintf(fid,'\n nobase = F');
        fprintf(fid,'\n dows=T');
        fprintf(fid,'\n WCONC = 1'); %Need to Gasparovic's correction later -- 15Apr2022 MGSaleh
        fprintf(fid,'\n sptype = ''mega-press-3''');
        fprintf(fid,'\n attmet = 1');
        fprintf(fid,'\n atth2o = 0.5'); %2');
        fprintf(fid,'\n PPMEND=1.7');
        fprintf(fid,'\n PPMST=4');
        fprintf(fid,'\n NAMREL = ''NAA+NAAG''');
        remove_met = {'PE', 'Asp', 'Cr','PCr','Asc','Ins','GPC','PCh','Scyllo','Tau','Ser','-CrCH2','Sing00','Sin215'};
        fprintf(fid,'\n NOMIT=%d', length(remove_met));
        for nmets = 1:length(remove_met)
            fprintf(fid,'\n CHOMIT(%d)=''%s''',nmets, remove_met{nmets});
        end
    case {'diff2'}
        if MRSCont.flags.isHERMES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERMES_Sig215_GSH.basis');
        elseif MRSCont.flags.isHERCULES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERC_Sig215_GSH.basis');
        end
        fprintf(fid,'\n atth2o = 0.5'); %4');
        fprintf(fid,'\n attmet = 1');
        fprintf(fid,'\n dows=T');
        fprintf(fid,'\n WCONC = 1'); %Need to Gasparovic's correction later -- 15Apr2022 MGSaleh
        fprintf(fid,'\n NAMREL = ''NAA+NAAG''');
        remove_met = {'PE', 'Asp', 'Cr','GABA','Ins','Gln','Glu','PCr','Scyllo','Tau','-CrCH2',...
                      'Ser','Asc','MM09','MM20','Lip09','Lip13a','Lip13b','Lip20','GPC','PCh','Sing00','Sing215'};
        fprintf(fid,'\n NOMIT=%d', length(remove_met));
        for nmets = 1:length(remove_met)
            fprintf(fid,'\n CHOMIT(%d)=''%s''',nmets, remove_met{nmets});
        end
        if MRSCont.flags.isHERMES==1
            fprintf(fid,'\n PPMEND=2'); %0.8 
            fprintf(fid,'\n PPMST=3.3'); %3.4
        elseif MRSCont.flags.isHERCULES==1
            fprintf(fid,'\n PPMEND=0.5');
            fprintf(fid,'\n PPMST=4');
        end
%         fprintf(fid,'\n WSMET = ''Sin215''');
%         fprintf(fid,'\n WSPPM = 2.15');
%         fprintf(fid,'\n N1HMET = 1');
    case {'sum'}
        if MRSCont.flags.isHERMES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERMES_SUM.basis');
        elseif MRSCont.flags.isHERCULES==1
            fprintf(fid,'\n FILBAS=''/Users/msaleh/Documents/MATLAB/LCModel/basis/%s''' ,'HERC_SUM.basis');
        end
        fprintf(fid,'\n attmet = 1');
        fprintf(fid,'\n atth2o = 0.25');
        fprintf(fid,'\n PPMEND=0.6');
        fprintf(fid,'\n PPMST=4');
        fprintf(fid,'\n dows=T');
        fprintf(fid,'\n WCONC = 1'); %Need to Gasparovic's correction later -- 15Apr2022 MGSaleh
        fprintf(fid,'\n NAMREL = ''Cr+PCr''');
        remove_met = {'Sing00','Sing215'};%,'Lip09'};
        fprintf(fid,'\n NOMIT=%d', length(remove_met));
        for nmets = 1:length(remove_met)
            fprintf(fid,'\n CHOMIT(%d)=''%s''',nmets, remove_met{nmets});
        end
    case {'A'}
        fprintf(fid,'\n attmet = 1');
        fprintf(fid,'\n atth2o = 1');
        fprintf(fid,'\n PPMEND=0.6');
        fprintf(fid,'\n PPMST=4');
        fprintf(fid,'\n dows=T');
        fprintf(fid,'\n WCONC = 1'); %Need to Gasparovic's correction later -- 15Apr2022 MGSaleh
        fprintf(fid,'\n NAMREL = ''Cr+PCr''');
        remove_met = {'Sing00'};%,'Lip09'};
        fprintf(fid,'\n NOMIT=%d', length(remove_met));
        for nmets = 1:length(remove_met)
            fprintf(fid,'\n CHOMIT(%d)=''%s''',nmets, remove_met{nmets});
        end
end
fprintf(fid,'\n $END');
fclose(fid);

% Save the paths to the control file and LCModel output files to MRSCont
controlWhich = strrep(subspec, 'outfile', 'controlfile');
outputWhich = strrep(subspec, 'outfile', 'outputfile');
MRSCont.opts.fit.lcmodel.(controlWhich){kk} = controlFile;
MRSCont.opts.fit.lcmodel.(outputWhich){kk} = outputFilesCommon;

end

function outputFilename = makeUniqueFileName(inputFilename)
    % For batch analysis, get the last two sub-folders (e.g. site and
    % subject) to augment the filename, avoiding duplicate output filenames
    [path,filename,ext]   = fileparts(inputFilename);
    path_split          = regexp(path,filesep,'split');
    
    if length(path_split) > 2
        outputFilename = [path_split{end-1} '_' path_split{end} '_' filename ext];
    else
        outputFilename = [path_split{end} '_' filename ext];
    end
end
