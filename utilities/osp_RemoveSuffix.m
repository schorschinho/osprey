function[OutName] = osp_RemoveSuffix(Name)
%% function[OutName] = osp_RemoveSuffix(Name)
% Function to remove suffix from filenames while respecting key-value pairs
% as defined in BIDs specification.
% i.e.  RemoveSuffix('MyDataSet_acq-press_act') = MyDataSet_acq-press
% while RemoveSuffix('MyDataSet_acq-press') = MyDataSet_acq-press

[~,Name,~] = fileparts(Name); %Remove lingering paths and extensions

SplitName = strsplit(Name,'_'); %Separate out by underscores

if isempty(strfind(SplitName{end},'-')) % If final underscore is not a key-value pair
    SplitName = SplitName(1:end-1); % Remove suffix
end
OutName = strjoin(SplitName,'_');

end