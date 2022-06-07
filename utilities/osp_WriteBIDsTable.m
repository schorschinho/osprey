function[] = osp_WriteBIDsTable(Table,OutLoc)
%% function[] = WriteBIDsTable(Table,OutLoc)
% Function for writing BIDs compliant tables.
%
%   USAGE:
%       WriteBIDsTable(Table,OutLoc,JSON)
%
%   INPUTS:
%       Table       = Matlab table structure, including json fields in
%                     table properties
%       OutLoc      = Path, including filename, but without extension

writetable(Table,[OutLoc,'.txt'],'Delimiter','\t'); % Write table with tab delimiter
movefile([OutLoc,'.txt'],[OutLoc,'.tsv']); % Change file extension to tsv

for JJ=1:length(Table.Properties.VariableNames)
    JSON.(Table.Properties.VariableNames{JJ}).LongName = Table.Properties.CustomProperties.VariableLongNames{JJ};
    JSON.(Table.Properties.VariableNames{JJ}).Description = Table.Properties.VariableDescriptions{JJ};
    JSON.(Table.Properties.VariableNames{JJ}).Units = Table.Properties.VariableUnits{JJ};
    %JSON.(Table.Properties.VariableNames{JJ}).TermURL = Table.Properties.CustomProperties.VariableTermURL{JJ};
end

% Write json structure to text file following jsonencode
fid=fopen([OutLoc,'.json'],'w');
if isMATLABReleaseOlderThan("R2021a")
    fprintf(fid, jsonencode(JSON)); 
else
    fprintf(fid, jsonencode(JSON,PrettyPrint=true)); %If newer Matlab version, utilize prettyprint formatting for JSON
end
fclose(fid);

end