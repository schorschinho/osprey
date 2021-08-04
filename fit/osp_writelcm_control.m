function MRSCont = osp_writelcm_control(MRSCont, kk, which, LCMparam)
%% [MRSCont] = osp_writelcm_control(MRSCont, kk, which, LCMparam)
%   This function creates a LCModel compatible .control file which can be used
%   for LCModel batch processing. 
%
%   USAGE:
%       RF = osp_writelcm_control(MRSCont,kk,which,LCMparam);
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
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
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
    LCMparam.DOWS = 'F';
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

% Create a control file with the same name as the water-suppressed file.
[~, rawFilename, ~] = fileparts(name_raw);
controlFile = fullfile(saveDestination, strrep(makeUniqueFileName(name_raw), '.RAW', '.control'));
% Create the LCModel output filenames (coord, print, ps, table, csv)
% Since the names are identical (and only the extension will vary), create
% the string without extension.
outputFiles = fullfile(MRSCont.outputFolder, 'LCMoutput', strrep(makeUniqueFileName(name_raw), '.RAW', ''));

fid = fopen(controlFile,'w+');
fprintf(fid,' $LCMODL');
fprintf(fid,'\n key = %i', LCMparam.key);
fprintf(fid,'\n OWNER = ''%s''', LCMparam.owner);
fprintf(fid,'\n Title = ''%s''', rawFilename);
fprintf(fid,'\n HZPPPM =%2.6e, DELTAT=%2.6e, NUNFIL=%i',MRSCont.processed.(which){kk}.txfrq/1e6, MRSCont.processed.(which){kk}.dwelltime, MRSCont.processed.(which){kk}.sz(1));
fprintf(fid,'\n FILBAS = ''%s''', MRSCont.opts.fit.basisSetFile);
fprintf(fid,'\n DKNTMN = %5.2f', LCMparam.DKNTMN);
fprintf(fid,'\n DOWS = %s', LCMparam.DOWS);
fprintf(fid,'\n FILRAW = ''%s''', name_raw);
if exist('name_ref','var')
    fprintf(fid,'\n FILH2O = ''%s''', selectedWaterRef);
    fprintf(fid,'\n ATTH2O = %5.2f', LCMparam.ATTH2O);
    fprintf(fid,'\n WCONC = %5.1f', LCMparam.WCONC);
end    
fprintf(fid,'\n ATTMET = %5.5f', LCMparam.ATTMET);
fprintf(fid,'\n NEACH = %i', LCMparam.NEACH);
fprintf(fid,'\n WDLINE(6) = %5.1f', LCMparam.WDLINE);
fprintf(fid,'\n PPMST = %5.2f', MRSCont.opts.fit.range(2));
fprintf(fid,'\n PPMEND = %5.2f', MRSCont.opts.fit.range(1));
fprintf(fid,'\n NSIMUL = %i', 12);
fprintf(fid,'\n NCOMBI = %i', length(LCMparam.CHCOMB));
for i = 1:length(LCMparam.CHCOMB)
    fprintf(fid,'\n CHCOMB(%i) = ''%s''', i, LCMparam.CHCOMB{i});
end
fprintf(fid,'\n NOMIT= %i',length(LCMparam.NOMIT));
for i = 1:length(LCMparam.NOMIT)
    fprintf(fid,'\n CHOMIT(%i) = ''%s''', i, LCMparam.NOMIT{i});
end
fprintf(fid,'\n NAMREL = ''%s''', LCMparam.NAMREL);
% fprintf(fid,'\n CONREL = %5.2f',LCMparam.CONREL);
fprintf(fid,'\n DOECC = %s', LCMparam.DOECC);
fprintf(fid,'\n LTABLE = %i', 7);
fprintf(fid,'\n LCSV = %i', 11);
fprintf(fid,'\n LCOORD = %i', 9);
fprintf(fid,'\n LPRINT = %i', 6);
fprintf(fid,'\n FILTAB = ''%s''', strcat(outputFiles, '.table'));
fprintf(fid,'\n FILCSV = ''%s''', strcat(outputFiles, '.csv'));
fprintf(fid,'\n FILCOO = ''%s''', strcat(outputFiles, '.coord'));
fprintf(fid,'\n FILPRI = ''%s''', strcat(outputFiles, '.print'));
fprintf(fid,'\n FILPS = ''%s''',  strcat(outputFiles, '.ps'));
fprintf(fid,'\n $END');
fclose(fid);

% Save the paths to the control file and LCModel output files to MRSCont
controlWhich = strrep(subspec, 'outfile', 'controlfile');
outputWhich = strrep(subspec, 'outfile', 'outputfile');
MRSCont.opts.fit.lcmodel.(controlWhich){kk} = controlFile;
MRSCont.opts.fit.lcmodel.(outputWhich){kk} = outputFiles;

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
