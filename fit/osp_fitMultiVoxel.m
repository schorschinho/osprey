function [outMRSCont] = osp_fitMultiVoxel(MRSCont)
metFitTime = tic;
outMRSCont= MRSCont;
fitMRSCont = MRSCont;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
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
                    fitMRSCont.fit.scale{kk} =  MRSCont.fit.scale{kk};
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
            outMRSCont.fit.resBasisSet = fitMRSCont.fit.resBasisSet;
            outMRSCont.fit.results = fitMRSCont.fit.results;
        else
            fields = {'resBasisSet','results'};
            for f = 1 : length(fields)
                if isfield(outMRSCont.fit,fields{f})
                    if iscell(outMRSCont.fit.(fields{f}))
                        %PRIAM data
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    else
                        temp = outMRSCont.fit.(fields{f});
                        outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
                        outMRSCont.fit.(fields{f}){1} = temp;
                        outMRSCont.fit.(fields{f}){x} = fitMRSCont.fit.(fields{f});
                    end
                end
            end
        end
    end
    time = toc(metFitTime);
    outMRSCont.runtime.FitMet = time;
elseif MRSCont.flags.isMRSI == 1
    vox = 1;
    %Fit center of MRSI first for inital guess parameters
    if isfield( MRSCont, 'mask')
        
        if max(max(max(MRSCont.mask{1}))) == 2
            [r, c] = find(MRSCont.mask{1}(:,:,2)== 2);
            cx = round(r);
            cy = round(c);
            cz =2;
            zList = [2,3,1];
        else
            [r, c] = find(MRSCont.mask{1}(:,:,2)== 1);
            cx = round(mean(c));
            cy = round(mean(r));
            cz =2;
            zList = [2,3,1];
        end
    else
       cx = round(XVox/2);
       cy = round(YVox/2);
       cz = 2;
    end

   for ss = 1 : NoSubSpec % Loop over Subspec
       for kk = 1 :MRSCont.nDatasets
           if ZVox <=1
               fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[cx,cy]);
           else
               fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[cx,cy,cz]);
           end
       end
   end
    fprintf('\nFitting center voxel (%d, %d, %d) of the MRSI dataset...\n', cx, cy, cz);

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
    outMRSCont.fit = fitMRSCont.fit;
    fields = {'resBasisSet','results'};
    for f = 1 : length(fields)
        if ZVox <=1
            temp = outMRSCont.fit.(fields{f});
            outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
            outMRSCont.fit.(fields{f}){cx,cy} = fitMRSCont.fit.(fields{f});
        else  % 3D MRSI data
            temp = outMRSCont.fit.(fields{f});
            outMRSCont.fit = rmfield(outMRSCont.fit, fields{f});
            outMRSCont.fit.(fields{f}){cx,cy,cz} = fitMRSCont.fit.(fields{f});
        end
    end

    for z = 1 : ZVox
        for x = 1 : XVox
            for y = 1 : YVox
                for f = 1 : length(fields)
                     % 2D MRSI data
                    if ZVox <=1
                        outMRSCont.fit.(fields{f}){x,y} = outMRSCont.fit.(fields{f}){cx,cy};
                    else  % 3D MRSI data
                        outMRSCont.fit.(fields{f}){x,y,z} = outMRSCont.fit.(fields{f}){cx,cy,cz};
                    end
                end
                if ~((x == cx) && (y == cy) && (z == cz)) % Do not overwrite the center voxel
                    if ZVox <=1
                        outMRSCont.fit.results{x,y}.off.fitParams{1}.ampl = zeros(size(outMRSCont.fit.results{x,y}.off.fitParams{1}.ampl));
                        outMRSCont.fit.results{x,y}.off.fitParams{1}.beta_j = zeros(size(outMRSCont.fit.results{x,y}.off.fitParams{1}.beta_j));
                        if MRSCont.flags.isMEGA
                            outMRSCont.fit.results{x,y}.diff1.fitParams{1}.ampl = zeros(size(outMRSCont.fit.results{x,y}.diff1.fitParams{1}.ampl));
                            outMRSCont.fit.results{x,y}.diff1.fitParams{1}.beta_j = zeros(size(outMRSCont.fit.results{x,y}.diff1.fitParams{1}.beta_j));
                        end
                    else
                        outMRSCont.fit.results{x,y,z}.off.fitParams{1}.ampl = zeros(size(outMRSCont.fit.results{x,y,z}.off.fitParams{1}.ampl));
                        outMRSCont.fit.results{x,y,z}.off.fitParams{1}.beta_j = zeros(size(outMRSCont.fit.results{x,y,z}.off.fitParams{1}.beta_j));
                        if MRSCont.flags.isMEGA
                            outMRSCont.fit.results{x,y,z}.diff1.fitParams{1}.ampl = zeros(size(outMRSCont.fit.results{x,y,z}.diff1.fitParams{1}.ampl));
                            outMRSCont.fit.results{x,y,z}.diff1.fitParams{1}.beta_j = zeros(size(outMRSCont.fit.results{x,y,z}.diff1.fitParams{1}.beta_j));
                        end
                    end
                end
            end
        end
    end

    % A.ph0(1) = fitMRSCont.fit.results.off.fitParams{kk}.ph0;
    % A.ph1(1)= fitMRSCont.fit.results.off.fitParams{kk}.ph1;
    % A.gaussLB(1) = fitMRSCont.fit.results.off.fitParams{kk}.gaussLB;
    % MRSCont.fit.MRSIpriors.A.ph0 = [];
    % MRSCont.fit.MRSIpriors.A.ph1 = [];
    % MRSCont.fit.MRSIpriors.A.gaussLB = [];
    % MRSCont.fit.MRSIpriors.A.prelimParams =fitMRSCont.fit.results.off.fitParams{kk}.prelimParams;
    % 
    % if MRSCont.flags.isMEGA
    %     diff1.ph0(1) = fitMRSCont.fit.results.diff1.fitParams{kk}.ph0;
    %     diff1.ph1(1)= fitMRSCont.fit.results.diff1.fitParams{kk}.ph1;
    %     diff1.gaussLB(1) = fitMRSCont.fit.results.diff1.fitParams{kk}.gaussLB;
    %     MRSCont.fit.MRSIpriors.diff1.ph0 = [];
    %     MRSCont.fit.MRSIpriors.diff1.ph1 = [];
    %     MRSCont.fit.MRSIpriors.diff1.gaussLB = [];
    %     MRSCont.fit.MRSIpriors.diff1.prelimParams =fitMRSCont.fit.results.diff1.fitParams{kk}.prelimParams;
    % end


    A.ph0(1) = fitMRSCont.fit.results.off.fitParams{kk}.ph0;
    A.ph1(1)= fitMRSCont.fit.results.off.fitParams{kk}.ph1;
    A.gaussLB(1) = fitMRSCont.fit.results.off.fitParams{kk}.gaussLB;
    MRSCont.fit.MRSIpriors.A.ph0 = A.ph0;
    MRSCont.fit.MRSIpriors.A.ph1 = A.ph1;
    MRSCont.fit.MRSIpriors.A.gaussLB = A.gaussLB;
    MRSCont.fit.MRSIpriors.A.prelimParams =fitMRSCont.fit.results.off.fitParams{kk}.prelimParams;

    if MRSCont.flags.isMEGA
        diff1.ph0(1) = fitMRSCont.fit.results.diff1.fitParams{kk}.ph0;
        diff1.ph1(1)= fitMRSCont.fit.results.diff1.fitParams{kk}.ph1;
        diff1.gaussLB(1) = fitMRSCont.fit.results.diff1.fitParams{kk}.gaussLB;
        MRSCont.fit.MRSIpriors.diff1.ph0 = diff1.ph0;
        MRSCont.fit.MRSIpriors.diff1.ph1 = diff1.ph1;
        MRSCont.fit.MRSIpriors.diff1.gaussLB = diff1.gaussLB;
        MRSCont.fit.MRSIpriors.diff1.prelimParams =fitMRSCont.fit.results.diff1.fitParams{kk}.prelimParams;
    end


    % Fit all remaining voxels
    for zi = 1 : ZVox
        z = zList(zi);
         for nVox = 1 : XVox*YVox
           [x,y] = osp_spiral(nVox);
           x = x+cx;
           y = y+cy;

            [~] = printLog('OspreyFit',[kk,vox],[MRSCont.nDatasets, XVox*YVox*ZVox],progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
            fprintf('\nFitting voxel (%d, %d, %d) of the MRSI dataset...\n', x, y, z); 
            try
             if MRSCont.mask{kk}(x,y,z)>0
               for ss = 1 : NoSubSpec % Loop over Subspec
                   for kk = 1 :MRSCont.nDatasets
                       if ZVox <=1
                           fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y]);
                       else
                           fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y,z]);
                       end
                   end
               end
                fitMRSCont.fit.MRSIpriors =  MRSCont.fit.MRSIpriors;  %Load prior results into the struct
                if ~((x == cx) && (y == cy) && (z == cz)) % Do not re-analyze the center voxel
                    fitMRSCont.fit.MRSIpriors.didFit = 0;
                    starttime = now; 
                    while (now<starttime+120/86400 ) && (fitMRSCont.fit.MRSIpriors.didFit == 0) 
                        if MRSCont.flags.isUnEdited
                            [fitMRSCont] = osp_fitUnEdited(fitMRSCont);
                            fitMRSCont.fit.MRSIpriors.didFit = 1;
                        elseif MRSCont.flags.isMEGA
                            [fitMRSCont] = osp_fitMEGA(fitMRSCont);
                            fitMRSCont.fit.MRSIpriors.didFit = 1;
                        elseif MRSCont.flags.isHERMES
                            [fitMRSCont] = osp_fitHERMES(fitMRSCont);
                            fitMRSCont.fit.MRSIpriors.didFit = 1;
                        elseif MRSCont.flags.isHERCULES
                            % For now, fit HERCULES like HERMES data
                            [fitMRSCont] = osp_fitHERCULES(fitMRSCont);
                            fitMRSCont.fit.MRSIpriors.didFit = 1;
                        else
                            msg = 'No flag set for sequence type!';
                            fprintf(fileID,msg);
                            error(msg);
                        end
                        fclose all
                        A.ph0(end+1) = fitMRSCont.fit.results.off.fitParams{kk}.ph0;
                        A.ph1(end+1)= fitMRSCont.fit.results.off.fitParams{kk}.ph1;
                        A.gaussLB(end+1) = fitMRSCont.fit.results.off.fitParams{kk}.gaussLB;
                        MRSCont.fit.MRSIpriors.A.ph0 = median(A.ph0);
                        MRSCont.fit.MRSIpriors.A.ph1=median(A.ph1);
                        MRSCont.fit.MRSIpriors.A.gaussLB = median(A.gaussLB);
                        % MRSCont.fit.MRSIpriors.A.ph0 = [];
                        % MRSCont.fit.MRSIpriors.A.ph1=[];
                        % MRSCont.fit.MRSIpriors.A.gaussLB = [];
                        if MRSCont.flags.isMEGA
                            diff1.ph0(end+1) = fitMRSCont.fit.results.diff1.fitParams{kk}.ph0;
                            diff1.ph1(end+1)= fitMRSCont.fit.results.diff1.fitParams{kk}.ph1;
                            diff1.gaussLB(end+1) = fitMRSCont.fit.results.diff1.fitParams{kk}.gaussLB;
                            MRSCont.fit.MRSIpriors.diff1.ph0 = median(diff1.ph0);
                            MRSCont.fit.MRSIpriors.diff1.ph1=median(diff1.ph1);
                            MRSCont.fit.MRSIpriors.diff1.gaussLB = median(diff1.gaussLB);
                            % MRSCont.fit.MRSIpriors.diff1.ph0 = [];
                            % MRSCont.fit.MRSIpriors.diff1.ph1=[];
                            % MRSCont.fit.MRSIpriors.diff1.gaussLB = [];
                        end
                        fields = {'resBasisSet','results'};
                        for f = 1 : length(fields)
                             % 2D MRSI data
                            if ZVox <=1
                                outMRSCont.fit.(fields{f}){x,y} = fitMRSCont.fit.(fields{f});
                            else  % 3D MRSI data
                                outMRSCont.fit.(fields{f}){x,y,z} = fitMRSCont.fit.(fields{f});
                            end
                        end
                        endtime = now;
                        outMRSCont.fit.time(x,y,z) = (endtime - starttime)*86400;
                    end
                end
             end
             catch
             end
             vox = vox+1;
        end
    end
    time = toc(metFitTime);
    outMRSCont.runtime.FitMet = time;
    outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
    outMRSCont.fit.scale = MRSCont.fit.scale;
end
    time = toc(metFitTime);
    [~] = printLog('MRSIdone',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
    outMRSCont.runtime.FitMet = time;
    outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
    outMRSCont.fit.scale = MRSCont.fit.scale;
end
