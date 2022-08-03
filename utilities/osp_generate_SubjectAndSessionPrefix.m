function[PreFix] = osp_generate_SubjectAndSessionPrefix(Name,kk)
%% function[PreFix] = osp_RemoveSuffix(Name)
% Function to add subject and session filenames while respecting key-value pairs
% as defined in BIDs specification.

try
    subs = num2str([contains(Name,{'sub'},'IgnoreCase',true) contains(Name,'S'+ digitsPattern,'IgnoreCase',true)]);
catch
    subs = num2str([contains(Name,{'sub'},'IgnoreCase',true) contains(Name,'S','IgnoreCase',true)]);
end

switch subs
    case '0  0'
        PreFix = sprintf('sub-%03d',kk);
    case '0  1'
        PreFix = sprintf('S-%03d',kk);
    otherwise
        PreFix = sprintf('sub-%03d',kk);
end

if ~contains(Name,{'ses'},'IgnoreCase',true)
    PreFix = [PreFix '_' sprintf('ses-%03d',1)];
else
    SepFiles =  split(Name, filesep);
    SepFiles(strcmp(SepFiles,''))=[];
    ind = find(contains(ses,'ses'));
    goUpNTimes = length(SepFiles) - ind;
    tempDir = fullfile(Name,repmat(['..' filesep],[1 goUpNTimes]));
    tempDir = dir(tempDir);
    tempDir = tempDir(1).folder;
    SepFiles =  split(tempDir, filesep);
    PreFix = [PreFix '_' SepFiles{end-1}];
end

end