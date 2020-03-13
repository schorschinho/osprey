function RF=osp_writelcm_control(MRSCont,kk,which,LCMparam);
%% [MRSCont] = osp_writelcm_control(MRSCont,kk,which,LCMparam)
%   This function creates a LCModel compatible .control file which can be used
%   for LCModel batch processing. 
%
%   USAGE:
%       RF=osp_writelcm_control(MRSCont,kk,which,LCMparam);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          =  Index for the kk-th dataset (optional. Default = 1)
%       which       = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%      LCMparam     = struct with the control file parameters this struct created
%                     by the LCMcontrol.m script included in the Osprey process folder
%                       
%
%   OUTPUTS:
%       RF          = control output.
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
%
%   HISTORY:
%       2020-02-11: First version of the code.
% Fall back to defaults if not provided
if nargin<4
    error('ERROR: no input LCModel control struct specified.  Aborting!!');    
    if nargin < 3
        which = 'off';
        if nargin < 2
            kk = 1;
            if nargin<1
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
        subspec = which;
    case {'diff1'}
        subspec = 'DIFF1';
    case {'diff2'}
        subspec = 'DIFF2';
    case {'sum'}
        subspec = 'SUM';
end
        

% For batch analysis, get the last two sub-folders (e.g. site and
% subject)
[path,filename,~]       = fileparts(MRSCont.files{kk});
path_split          = regexp(path,filesep,'split');
if length(path_split) > 2
    name = [path_split{end-1} '_' path_split{end} '_' filename];
end
nameraw = [name '_LCM_' subspec];
% For batch analysis, get the last two sub-folders (e.g. site and
% subject)
if strcmpi(MRSCont.vendor, 'GE')
    [path_ref,filename_ref,~]   = fileparts(MRSCont.files{kk});
else
    [path_ref,filename_ref,~]   = fileparts(MRSCont.files_ref{kk});
end
path_ref_split          = regexp(path_ref,filesep,'split');
if length(path_ref_split) > 2
    name_ref = [path_ref_split{end-1} '_' path_ref_split{end} '_' filename_ref];
end


%write to txt file
fid=fopen(fullfile(saveDestination,[nameraw '.control']),'w+');
fprintf(fid,' $LCMODL');
fprintf(fid,'\n key = %i', LCMparam.key);
fprintf(fid,'\n OWNER = ''%s''', LCMparam.owner);
fprintf(fid,'\n Title = ''%s''', filename);
fprintf(fid,'\n HZPPPM =%2.6e, DELTAT=%2.6e, NUNFIL=%i',MRSCont.processed.(which){kk}.txfrq/1e6,MRSCont.processed.(which){kk}.dwelltime,MRSCont.processed.(which){kk}.sz(1));
fprintf(fid,'\n FILBAS = ''%s''', LCMparam.FILBAS);
fprintf(fid,'\n DKNTMN = %5.2f',LCMparam.DKNTMN);
fprintf(fid,'\n DOWS = %s',LCMparam.DOWS);
fprintf(fid,'\n FILH2O = ''%s''', [LCMparam.FOLDER '/ref/' name_ref '_LCM_REF.RAW']);
fprintf(fid,'\n ATTH2O = %5.2f',LCMparam.ATTH2O);
fprintf(fid,'\n ATTMET = %5.5f',LCMparam.ATTMET);
fprintf(fid,'\n WCONC = %5.1f',LCMparam.WCONC);
fprintf(fid,'\n NEACH = %i',LCMparam.NEACH);
fprintf(fid,'\n WDLINE(6) = %5.1f',LCMparam.WDLINE);
fprintf(fid,'\n PPMST = %5.2f',MRSCont.opts.fit.range(2));
fprintf(fid,'\n PPMEND = %5.2f',MRSCont.opts.fit.range(1));
fprintf(fid,'\n NSIMUL = %i',12);
fprintf(fid,'\n NCOMBI = %i',length(LCMparam.CHCOMB));
for i = 1 : length(LCMparam.CHCOMB)
    fprintf(fid,'\n CHCOMB(%i) = ''%s''',i, LCMparam.CHCOMB{i});
end
fprintf(fid,'\n NOMIT= %i',length(LCMparam.NOMIT));
for i = 1 : length(LCMparam.NOMIT)
    fprintf(fid,'\n CHOMIT(%i) = ''%s''',i, LCMparam.NOMIT{i});
end
fprintf(fid,'\n NAMREL = ''%s''',LCMparam.NAMREL);
fprintf(fid,'\n CONREL = %5.2f',LCMparam.CONREL);
fprintf(fid,'\n DOECC = %s',LCMparam.DOECC);
fprintf(fid,'\n LTABLE = %i',7);
fprintf(fid,'\n LCSV = %i',11);
fprintf(fid,'\n LCOORD = %i',9);
fprintf(fid,'\n LPRINT = %i',6);
fprintf(fid,'\n FILTAB = ''%s''',[LCMparam.FOLDER '/LCMoutput/' nameraw '.table']);
fprintf(fid,'\n FILCSV = ''%s''',[LCMparam.FOLDER '/LCMoutput/' nameraw '.csv']);
fprintf(fid,'\n FILCOO = ''%s''',[LCMparam.FOLDER '/LCMoutput/' nameraw '.coord']);
fprintf(fid,'\n FILPRI = ''%s''',[LCMparam.FOLDER '/LCMoutput/' nameraw '.print']);
fprintf(fid,'\n FILRAW = ''%s''',[LCMparam.FOLDER '/metabs/' nameraw '.RAW']);
fprintf(fid,'\n FILPS = ''%s''',[LCMparam.FOLDER '/LCMoutput/' nameraw '.PS']);
fprintf(fid,'\n $END');
fclose(fid);
