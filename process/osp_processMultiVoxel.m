function [outMRSCont] = osp_processMultiVoxel(MRSCont)    
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
outMRSCont= MRSCont;
procMRSCont = MRSCont;
%% Get infos to set up a loop to process all voxels
if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
else if MRSCont.flags.isMRSI == 1
        XVox = MRSCont.raw{1}.nXvoxels;
        YVox = MRSCont.raw{1}.nYvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end


if MRSCont.flags.isPRIAM == 1
for x = 1 : XVox
    for kk = 1 :MRSCont.nDatasets
        procMRSCont.raw{kk} = op_takeVoxel(MRSCont.raw{kk},x);
        if procMRSCont.flags.hasRef
            procMRSCont.raw_ref{kk} = op_takeVoxel(MRSCont.raw_ref{kk},x);
        end
        if procMRSCont.flags.hasWater
            procMRSCont.raw_w{kk} = op_takeVoxel(MRSCont.raw_w{kk},x);
        end
    end
    if MRSCont.flags.isUnEdited
        [procMRSCont] = osp_processUnEdited(procMRSCont);
    elseif MRSCont.flags.isMEGA           
        [procMRSCont] = osp_processMEGA(procMRSCont);
    elseif MRSCont.flags.isHERMES
        [procMRSCont] = osp_processHERMES(procMRSCont);
    elseif MRSCont.flags.isHERCULES
        % For now, process HERCULES like HERMES data
        [procMRSCont] = osp_processHERCULES(procMRSCont);
    else
        msg = 'No flag set for sequence type!';
        fprintf(fileID,msg);
        error(msg);
    end  
    if x == 1
        outMRSCont.processed = procMRSCont.processed;
        outMRSCont.QM = procMRSCont.QM;
    else
       SubSpecNames = fieldnames(outMRSCont.processed);
       NoSubSpec = length(fieldnames(outMRSCont.processed)); 
       for ss = 1 : NoSubSpec % Loop over Subspec 
           for kk = 1 :MRSCont.nDatasets
                outMRSCont.processed.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed.(SubSpecNames{ss}){kk},procMRSCont.processed.(SubSpecNames{ss}){kk},x);                
           end
       end
       %Adding QM fields
        fields = {'QM'};
        for f = 1 : length(fields)
            if isfield(outMRSCont,fields{f})
                if iscell(outMRSCont.(fields{f}))
                    %PRIAM data
%                             if length(index)==1
                        outMRSCont.(fields{f}){x} = procMRSCont.(fields{f});
%                             end

%                             % 2D MRSI data
%                             if length(index)==2
%                                 outMRSCont.(fields{f}){index(1),index(2)} = procMRSCont.(fields{f});
%                             end
% 
%                             % 3D MRSI data
%                             if length(index)==3
%                                 outMRSCont.(fields{f}){index(1),index(2),index(3)} = procMRSCont.(fields{f});
%                             end
                else
%                             if length(index)==1
                        temp = outMRSCont.(fields{f});
                        outMRSCont = rmfield(outMRSCont, fields{f});
                        outMRSCont.(fields{f}){1} = temp;
                        outMRSCont.(fields{f}){x} = procMRSCont.(fields{f});
%                             end

                    % 2D MRSI data
%                             if length(index)==2
%                                 temp = outMRSCont.(fields{f});
%                                 out = rmfield(out, fields{f});
%                                 outMRSCont.(fields{f}){1} = temp;
%                                 outMRSCont.(fields{f}){index(1),index(2)} = procMRSCont.(fields{f});
%                             end

                    % 3D MRSI data
%                             if length(index)==3
%                                 temp = outMRSCont.(fields{f});
%                                 out = rmfield(out, fields{f});
%                                 outMRSCont.(fields{f}){1} = temp;
%                                 outMRSCont.(fields{f}){index(1),index(2),index(3)} = procMRSCont.(fields{f});
%                             end
                end            
            end
        end
    end
end
    
    
end