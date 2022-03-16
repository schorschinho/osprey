function[] = osp_WriteBIDsTable(Table,OutLoc,JSON)
%% function[] = WriteBIDsTable(Table,OutLoc,JSON)
% Function for writing BIDs compliant tables. Optional variable, JSON for
% writing the sidecar file.
%
%   USAGE:
%       WriteBIDsTable(Table,OutLoc,JSON)
%
%   INPUTS:
%       Table       = Matlab table structure
%       OutLoc      = Path, including filename, but without extension
%       JSON        = (Optional) Structure to be encoded in JSON sidecar

writetable(Table,[OutLoc,'.txt'],'Delimiter','\t'); % Write table with tab delimiter
movefile([OutLoc,'.txt'],[OutLoc,'.tsv']); % Change file extension to tsv

% If JSON structure included, write this to a sidecar file with the same
% name and location, but with .json extension.
if exist('JSON','var')
    fid=fopen([OutLoc,'.json'],'w');
    fprintf(fid, jsonencode(JSON,'PrettyPrint',true)); 
    fclose(fid);
end

end