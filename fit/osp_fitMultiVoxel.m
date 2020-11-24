function [outMRSCont] = osp_fitMultiVoxel(MRSCont)    
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
metFitTime = tic;
outMRSCont= MRSCont;
fitMRSCont = MRSCont;
%% Get infos to set up a loop to process all voxels
if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
else if MRSCont.flags.isMRSI == 1
        XVox = MRSCont.raw{1}.nXvoxels;
        YVox = MRSCont.raw{1}.nYvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end
SubSpecNames = fieldnames(fitMRSCont.processed);
NoSubSpec = length(fieldnames(fitMRSCont.processed)); 

if MRSCont.flags.isPRIAM == 1
for x = 1 : XVox
   for ss = 1 : NoSubSpec % Loop over Subspec 
        for kk = 1 :MRSCont.nDatasets
                fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},x);                
       end
   end
    for kk = 1 :MRSCont.nDatasets % Loop over scale values
                fitMRSCont.fit.scale{kk} =  MRSCont.fit.scale{kk}(x);                
    end
    if MRSCont.flags.isUnEdited
        [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
    elseif MRSCont.flags.isMEGA
        [fitMRSCont] = osp_fitMEGA(fitMRSCont);
    elseif MRSCont.flags.isHERMES
        [fitMRSCont] = osp_fitHERMES(fitMRSCont);
    elseif MRSCont.flags.isHERCULES
        % For now, fit HERCULES like HERMES data
        [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
    else
        msg = 'No flag set for sequence type!';
        fprintf(fileID,msg);
        error(msg);
    end
    
    if x == 1
        outMRSCont.fit = fitMRSCont.fit;
    else
        fields = {'resBasisSet','results'};
        for f = 1 : length(fields)
            if isfield(outMRSCont.fit,fields{f})
                if iscell(outMRSCont.fit.(fields{f}))
                    %PRIAM data
%                             if length(index)==1
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
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
                        temp = outMRSCont.fit.(fields{f});
                        outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
                        outMRSCont.fit.(fields{f}){1} = temp;
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
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
time = toc(metFitTime);
outMRSCont.runtime.FitMet = time;    
outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
outMRSCont.fit.scale = MRSCont.fit.scale;
end