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
fprintf('Writing table to file = %s\n', [OutLoc,'.tsv']);

if length(Table.Properties.VariableNames) > 0
	for JJ=1:length(Table.Properties.VariableNames)
		JSON.(Table.Properties.VariableNames{JJ}).LongName = Table.Properties.CustomProperties.VariableLongNames{JJ};
		JSON.(Table.Properties.VariableNames{JJ}).Description = Table.Properties.VariableDescriptions{JJ};
		JSON.(Table.Properties.VariableNames{JJ}).Units = Table.Properties.VariableUnits{JJ};
		%JSON.(Table.Properties.VariableNames{JJ}).TermURL = Table.Properties.CustomProperties.VariableTermURL{JJ};
	end
	
	% Write json structure to text file following jsonencode
	json_fn = [OutLoc,'.json'];
	fprintf('Writing table to file = %s\n', json_fn);
	fid=fopen(json_fn,'w');
	fprintf(fid, jsonencode(JSON)); 
	fclose(fid);
end

end