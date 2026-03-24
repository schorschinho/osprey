function [outMRSCont] = osp_fitMultiVoxel_IDEAL_Julia(MRSCont)
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
        XVox = MRSCont.raw{1}.nYvoxels;
        YVox = MRSCont.raw{1}.nXvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end
SubSpecNames = fieldnames(fitMRSCont.processed);
NoSubSpec = length(fieldnames(fitMRSCont.processed));

steps = [20 24;
        10 12;
        5 6;
        2 3;
        1 1;];
    

dims_steps = [20 24;
        2 2;
        4 4;
        10 12;
        20 24];    
for kk = 1 : MRSCont.nDatasets
%Outer loop with resampling steps
for res = 1 : 5
    %Resample brain mask
    mask = imresize(squeeze(double(fitMRSCont.mask{kk}(:,:,2))),[dims_steps(res,2) dims_steps(res,1)]);
    fitMRSCont = MRSCont;
    for xx = 1 : fitMRSCont.processed.A{kk}.nXvoxels
        for yy = 1 : fitMRSCont.processed.A{kk}.nYvoxels
            for zz = 1 : fitMRSCont.processed.A{kk}.nZvoxels
                if fitMRSCont.mask{kk}(yy,xx,zz) == 0
                    fitMRSCont.processed.A{kk}.specs(:,xx,yy,zz) = ones(512,1) * nan;
                    fitMRSCont.processed.A{kk}.fids(:,xx,yy,zz) = ones(512,1) * nan;
                end
            end
        end
    end
    if res == 1                
        spec = nanmean(fitMRSCont.processed.A{1, 1}.specs,2);
        spec = nanmean(spec,3);
        spec = squeeze(nanmean(spec,4));
        
        fitMRSCont.processed.A{1, 1}.nXvoxels=1;
        fitMRSCont.processed.A{1, 1}.nYvoxels=1;
        fitMRSCont.processed.A{1, 1}.nZvoxels=1;
        fitMRSCont.processed.A{1,1}.specs = spec;
        fitMRSCont.processed.A{1,1}.fids = ifft(fftshift(fitMRSCont.processed.A{1,1}.specs,1),[],1);

        [fitMRSCont] = osp_fitUnEdited(fitMRSCont);


%         for z = 1 : ZVox
%             for x = 1 : XVox
%                 for y = 1 : YVox
%                     for f = 1 : length(fields)
%                          % 2D MRSI data
%                         if ZVox <=1
%                             outMRSCont.fit.(fields{f}){x,y} = outMRSCont.fit.(fields{f}){cx,cy};
%                         else  % 3D MRSI data
%                             outMRSCont.fit.(fields{f}){x,y,z} = outMRSCont.fit.(fields{f}){cx,cy,cz};
%                         end
%                     end
%                     if ~((x == cx) && (y == cy) && (z == cz)) % Do not overwrite the center voxel
%                         outMRSCont.fit.results{x,y}.off.fitParams{1}.ampl = zeros(size(outMRSCont.fit.results{x,y}.off.fitParams{1}.ampl));
%                         outMRSCont.fit.results{x,y}.off.fitParams{1}.beta_j = zeros(size(outMRSCont.fit.results{x,y}.off.fitParams{1}.beta_j));
%                     end
%                 end
%             end
%         end
        ph0(1) = fitMRSCont.fit.results.off.fitParams{kk}.ph0;
        ph1(1)= fitMRSCont.fit.results.off.fitParams{kk}.ph1;
        gaussLB(1) = fitMRSCont.fit.results.off.fitParams{kk}.gaussLB;

       MRSCont.fit.MRSIfitPriors = fitMRSCont.fit;
    else
    % Create mean data with new dimensions 


    temp_spec = mat2cell(fitMRSCont.processed.A{1, 1}.specs,512,[steps(1,res) steps(1,res)], [steps(2,res) steps(2,res)]);
    for newXVox = 1 : dims_steps(res,1)
        for newYVox = 1 : dims_steps(res,2)
            spec(:,newXVox,newYVox) =squeeze(mean(mean(cell2mat(temp_spec(1,newXVox,newYVox)),2),3));
        end
    end

    fitMRSCont.processed.A{1, 1}.nXvoxels=dims_steps(1,res);
    fitMRSCont.processed.A{1, 1}.nYvoxels=dims_steps(2,res);
    fitMRSCont.processed.A{1, 1}.nZvoxels=1;
    fitMRSCont.processed.A{1,1}.specs = spec;
    fitMRSCont.processed.A{1,1}.fids = ifft(fftshift(fitMRSCont.processed.A{1,1}.specs,1),[],1);
    for y = 1 : YVox
         for x = 1 : XVox
             try
             if mask(x,y)
               for ss = 1 : NoSubSpec % Loop over Subspec
                   for kk = 1 :MRSCont.nDatasets
                       if ZVox <=1
                           fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y]);
                       else
                           fitMRSCont.processed.(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},[x,y,z]);
                       end
                   end
               end
                fitMRSCont.fit =  MRSCont.fit.MRSIfitPriors;  %Load prior results into the struct
                [fitMRSCont] = osp_fitUnEdited(fitMRSCont);

                    fclose all
                    ph0(end+1) = fitMRSCont.fit.results.off.fitParams{kk}.ph0;
                    ph1(end+1)= fitMRSCont.fit.results.off.fitParams{kk}.ph1;
                    gaussLB(end+1) = fitMRSCont.fit.results.off.fitParams{kk}.gaussLB;
                    MRSCont.fit.MRSIfitPriors = fitMRSCont.fit; %Store new priors
                    MRSCont.fit.MRSIfitPriors.results.off.fitParams{kk}.ph0 = median(ph0);
                    MRSCont.fit.MRSIfitPriors.results.off.fitParams{kk}.ph1=median(ph1);
                    MRSCont.fit.MRSIfitPriors.results.off.fitParams{kk}.gaussLB = median(gaussLB);

                end
             
             catch
                 end
         end
        end
end
end
    time = toc(metFitTime);
    outMRSCont.runtime.FitMet = time;
    outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
    outMRSCont.fit.scale = MRSCont.fit.scale;
    
    
    time = toc(metFitTime);
    [~] = printLog('MRSIdone',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
    outMRSCont.runtime.FitMet = time;
    outMRSCont.fit.basisSet = MRSCont.fit.basisSet;
    outMRSCont.fit.scale = MRSCont.fit.scale;
end
