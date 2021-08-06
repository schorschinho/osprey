function MRSCont = osp_writelcm_control(MRSCont, kk, rr, which, LCMparam)
%% [MRSCont] = osp_writelcm_control(MRSCont, kk, rr, which, LCMparam)
%   This function creates a LCModel compatible .control file which can be used
%   for LCModel batch processing. 
%
%   USAGE:
%       RF = osp_writelcm_control(MRSCont, kk, rr, which, LCMparam);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = Index for the kk-th dataset (optional. Default = 1)
%       rr          = Index for the rr-th sub-spectrum of the kk-th dataset
%                     (optional. Default = 1)
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
if nargin < 5
    error('ERROR: no input LCModel control struct specified.  Aborting!!');
    if nargin < 4
        which = 'A';
        if nargin < 3
            rr = 1;
            if nargin < 2
                kk = 1;
                if nargin < 1
                    error('ERROR: no input Osprey container specified.  Aborting!!');
                end
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
name_raw = MRSCont.opts.fit.lcmodel.(subspec){kk}{rr};

% Retrieve the filenames of the LCModel .RAW files that were created at the
% end of OspreyProcess using osp_saveLCM.

if isfield(MRSCont.opts.fit.lcmodel, 'outfileRef')
    name_ref    = MRSCont.opts.fit.lcmodel.outfileRef{kk}{rr};
end
if isfield(MRSCont.opts.fit.lcmodel, 'outfileW')
    name_w      = MRSCont.opts.fit.lcmodel.outfileW{kk}{rr};
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
outputFilesCommon = fullfile(MRSCont.outputFolder, 'LCMoutput', strrep(makeUniqueFileName(name_raw), '.RAW', ''));

% Add individual LCM control parameters:
LCMparam = osp_editControlParameters(LCMparam, 'filraw', ['''' name_raw '''']);     % path to the .RAW file produced previously
LCMparam = osp_editControlParameters(LCMparam, 'filps',  ['''' strcat(outputFilesCommon, '.ps') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filtab', ['''' strcat(outputFilesCommon, '.table') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filcsv', ['''' strcat(outputFilesCommon, '.csv') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filcoo', ['''' strcat(outputFilesCommon, '.coord') '''']);
LCMparam = osp_editControlParameters(LCMparam, 'filpri', ['''' strcat(outputFilesCommon, '.print') '''']);

if exist('name_ref','var')
    LCMparam = osp_editControlParameters(LCMparam, 'filh2o',  ['''' selectedWaterRef '''']);
end    

% Write every line into the control file
fid = fopen(controlFile,'w+');
fprintf(fid,' $LCMODL');
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
fprintf(fid,'\n $END');
fclose(fid);

% Save the paths to the control file and LCModel output files to MRSCont
controlWhich = strrep(subspec, 'outfile', 'controlfile');
outputWhich = strrep(subspec, 'outfile', 'outputfile');
MRSCont.opts.fit.lcmodel.(controlWhich){kk}{rr} = controlFile;
MRSCont.opts.fit.lcmodel.(outputWhich){kk}{rr} = outputFilesCommon;

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
