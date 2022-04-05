function [outMRSCont] = osp_processMultiVoxel(MRSCont)    
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
outMRSCont= MRSCont;
procMRSCont = MRSCont;
if ~strcmp(MRSCont.opts.MoCo, 'none')
    procMRSCont_no_MoCo =MRSCont;
end
%% Get infos to set up a loop to process all voxels
time = tic;
if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
else if MRSCont.flags.isMRSI == 1
        XVox = MRSCont.raw{1}.nXvoxels;
        YVox = MRSCont.raw{1}.nYvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end

if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
if MRSCont.flags.isPRIAM == 1
    for x = 1 : XVox
        for kk = 1 :MRSCont.nDatasets
             [~] = printLog('OspreyProcess',[kk,x],[MRSCont.nDatasets, XVox],progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
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
            outMRSCont.runtime.Proc = procMRSCont.runtime.Proc;
            outMRSCont.QM = procMRSCont.QM;
        else
           SubSpecNames = fieldnames(outMRSCont.processed);
           NoSubSpec = length(fieldnames(outMRSCont.processed)); 
           for ss = 1 : NoSubSpec % Loop over Subspec 
               for kk = 1 :MRSCont.nDatasets
                    outMRSCont.processed.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed.(SubSpecNames{ss}){kk},procMRSCont.processed.(SubSpecNames{ss}){kk},x);                
               end
           end
           outMRSCont.runtime.Proc = outMRSCont.runtime.Proc + procMRSCont.runtime.Proc;
           %Adding QM fields
            fields = {'QM'};
            for f = 1 : length(fields)
                if isfield(outMRSCont,fields{f})
                    if iscell(outMRSCont.(fields{f}))
                            outMRSCont.(fields{f}){x} = procMRSCont.(fields{f});
                    else
                            temp = outMRSCont.(fields{f});
                            outMRSCont = rmfield(outMRSCont, fields{f});
                            outMRSCont.(fields{f}){1} = temp;
                            outMRSCont.(fields{f}){x} = procMRSCont.(fields{f});
                    end            
                end
            end
        end
    end    

elseif MRSCont.flags.isMRSI == 1    
    NVox = XVox*YVox*ZVox;
    vox = 1;
    for x = 1 : XVox
        for y = 1 : YVox
            for z = 1 : ZVox 
                for kk = 1 :MRSCont.nDatasets
                    [~] = printLog('OspreyProcess',[kk,vox],[MRSCont.nDatasets, NVox],progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
                    if ZVox <=1
                        procMRSCont.raw{kk} = op_takeVoxel(MRSCont.raw{kk},[x y]);
                        if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                            procMRSCont_no_MoCo.raw{kk} = op_takeVoxel(MRSCont.raw_no_MoCo{kk},[x y]);
                        end
                    else
                        procMRSCont.raw{kk} = op_takeVoxel(MRSCont.raw{kk},[x y z]);
                        if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                            procMRSCont_no_MoCo.raw{kk} = op_takeVoxel(MRSCont.raw_no_MoCo{kk},[x y z]);
                        end
                    end
                    if procMRSCont.flags.hasRef
                        try
                            if ZVox <=1
                                procMRSCont.raw_ref{kk} = op_takeVoxel(MRSCont.raw_ref{kk},[x y]);
                                if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                    procMRSCont_no_MoCo.raw_ref{kk} = op_takeVoxel(MRSCont.raw_ref{kk},[x y]);
                                end
                            else
                                 procMRSCont.raw_ref{kk} = op_takeVoxel(MRSCont.raw_ref{kk},[x y z]);
                                 if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                     procMRSCont_no_MoCo.raw_ref{kk} = op_takeVoxel(MRSCont.raw_ref{kk},[x y z]);
                                 end
                            end
                        catch
                        end
                    end
                    if procMRSCont.flags.hasWater
                        try
                            if ZVox <=1
                                procMRSCont.raw_w{kk} = op_takeVoxel(MRSCont.raw_w{kk},[x y]);
                                if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                    procMRSCont_no_MoCo.raw_w{kk} = op_takeVoxel(MRSCont.raw_w{kk},[x y]);
                                end
                            else
                                procMRSCont.raw_w{kk} = op_takeVoxel(MRSCont.raw_w{kk},[x y z]);
                                if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                     procMRSCont_no_MoCo.raw_w{kk} = op_takeVoxel(MRSCont.raw_w{kk},[x y z]);
                                end
                            end
                        catch
                        end
                    end
                end
                if MRSCont.flags.isUnEdited
                    [procMRSCont] = osp_processUnEdited(procMRSCont);
                    if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                        [procMRSCont_no_MoCo] = osp_processUnEdited(procMRSCont_no_MoCo);
                    end
                elseif MRSCont.flags.isMEGA           
                    [procMRSCont] = osp_processMEGA(procMRSCont);
                    if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                        [procMRSCont_no_MoCo] = osp_processMEGA(procMRSCont_no_MoCo);
                    end
                elseif MRSCont.flags.isHERMES
                    [procMRSCont] = osp_processHERMES(procMRSCont);
                    if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                        [procMRSCont_no_MoCo] = osp_processHERMES(procMRSCont_no_MoCo);
                    end
                elseif MRSCont.flags.isHERCULES
                    % For now, process HERCULES like HERMES data
                    [procMRSCont] = osp_processHERCULES(procMRSCont);
                    if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                        [procMRSCont_no_MoCo] = osp_processHERCULES(procMRSCont_no_MoCo);
                    end
                else
                    msg = 'No flag set for sequence type!';
                    fprintf(fileID,msg);
                    error(msg);
                end  
                if (x == 1) && (y == 1) && (z == 1)
                    outMRSCont.processed = procMRSCont.processed;
                    outMRSCont.QM = procMRSCont.QM;
                    if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                        outMRSCont.processed_no_MoCo = procMRSCont_no_MoCo.processed;
                        outMRSCont.QM_no_MoCo = procMRSCont_no_MoCo.QM;
                    end
                else
                    SubSpecNames = fieldnames(outMRSCont.processed);
                    NoSubSpec = length(fieldnames(outMRSCont.processed)); 
                    for ss = 1 : NoSubSpec % Loop over Subspec 
                        for kk = 1 :MRSCont.nDatasets
                            if ZVox <=1
                                outMRSCont.processed.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed.(SubSpecNames{ss}){kk},procMRSCont.processed.(SubSpecNames{ss}){kk},[x y]);
                            else
                                outMRSCont.processed.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed.(SubSpecNames{ss}){kk},procMRSCont.processed.(SubSpecNames{ss}){kk},[x y z]);
                            end
                            if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                if ZVox <=1
                                    outMRSCont.processed_no_MoCo.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed_no_MoCo.(SubSpecNames{ss}){kk},procMRSCont_no_MoCo.processed.(SubSpecNames{ss}){kk},[x y]);
                                else
                                    outMRSCont.processed_no_MoCo.(SubSpecNames{ss}){kk} = op_addVoxel(outMRSCont.processed_no_MoCo.(SubSpecNames{ss}){kk},procMRSCont_no_MoCo.processed.(SubSpecNames{ss}){kk},[x y z]);
                                end
                            end
                        end
                    end
                        %Adding QM fields
                        fields = {'QM'};
                        for f = 1 : length(fields)
                            if isfield(outMRSCont,fields{f})
                                if iscell(outMRSCont.(fields{f}))
                                    % 2D MRSI data
                                    if ZVox <=1
                                        outMRSCont.(fields{f}){x,y} = procMRSCont.(fields{f});
                                    else % 3D MRSI data
                                        outMRSCont.(fields{f}){x,y,z} = procMRSCont.(fields{f});
                                    end
                                else
                                    % 2D MRSI data
                                    if ZVox <=1
                                        temp = outMRSCont.(fields{f});
                                        outMRSCont = rmfield(outMRSCont, fields{f});
                                        outMRSCont.(fields{f}){1,1} = temp;
                                        outMRSCont.(fields{f}){x,y} = procMRSCont.(fields{f});
                                    else % 3D MRSI data
                                        temp = outMRSCont.(fields{f});
                                        outMRSCont = rmfield(outMRSCont, fields{f});
                                        outMRSCont.(fields{f}){1,1,1} = temp;
                                        outMRSCont.(fields{f}){x,y,z} = procMRSCont.(fields{f});
                                    end
                                end %Initial cell aray set up           
                            end %inital struct set up
                            if ~strcmp(MRSCont.opts.MoCo.target, 'none')
                                if isfield(outMRSCont,fields{f})
                                    if iscell(outMRSCont.([fields{f} '_no_MoCo']))
                                        % 2D MRSI data
                                        if ZVox <=1
                                            outMRSCont.([fields{f} '_no_MoCo']){x,y} = procMRSCont_no_MoCo.(fields{f});
                                        else % 3D MRSI data
                                            outMRSCont.([fields{f} '_no_MoCo']){x,y,z} = procMRSCont_no_MoCo.(fields{f});
                                        end
                                    else
                                        % 2D MRSI data
                                        if ZVox <=1
                                            temp = outMRSCont.(fields{f});
                                            outMRSCont = rmfield(outMRSCont, [fields{f} '_no_MoCo']);
                                            outMRSCont.([fields{f} '_no_MoCo']){1,1} = temp;
                                            outMRSCont.([fields{f} '_no_MoCo']){x,y} = procMRSCont_no_MoCo.(fields{f});
                                        else % 3D MRSI data
                                            temp = outMRSCont.(fields{f});
                                            outMRSCont = rmfield(outMRSCont, [fields{f} '_no_MoCo']);
                                            outMRSCont.([fields{f} '_no_MoCo']){1,1,1} = temp;
                                            outMRSCont.([fields{f} '_no_MoCo']){x,y,z} = procMRSCont_no_MoCo.(fields{f});
                                        end
                                    end %Initial cell aray set up           
                                end %inital struct set up
                            end
                        end % Fields in QM struct
                end % Initial set up for x=y=z=1
                vox = vox +1;
            end %zVox
        end %yVox
    end % xVox
end %isMRSI  
[~] = printLog('MRSIdone',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
end