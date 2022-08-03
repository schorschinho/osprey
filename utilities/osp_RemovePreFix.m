function[OutName] = osp_RemovePreFix(Name)
%% function[OutName] = osp_RemovePreFix(Name)
% Function to remove prefix from filenames while respecting key-value pairs
% as defined in BIDs specification.
% i.e.  osp_RemovePreFix('sub-001_MyDataSet_acq-press_act') = MyDataSet_acq-press
% while osp_RemovePreFix('sub-001_MyDataSet_acq-press') = sub-001_MyDataSet_acq-press

[~,Name,~] = fileparts(Name); %Remove lingering paths and extensions

SplitName = strsplit(Name,'_'); %Separate out by underscores

% if isempty(strfind(SplitName{1},'-')) % If final underscore is not a key-value pair
    SplitName = SplitName(2:end); % Remove Pre
% end
OutName = strjoin(SplitName,'_');

end