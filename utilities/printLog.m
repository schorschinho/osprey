function [msg] = printLog(Module,kk,nDatasets,progressText,GUI,MRSI)
%% [MRSCont] = printLog(MRSCont)
%   This function allows you to rebase your derviatives folder and files in
%   case you have processed them on another machine
%
%   USAGE:
%       MRSCont = OspreyRebase(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2021-05-06)
%       hzoelln2@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-05-06: First version of the code.


msg = '';
%% Load
 
if strcmp(Module,'OspreyLoad')
     if kk == 1    
        msg = sprintf('Loading raw data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
     msg = sprintf('Loading raw data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, nDatasets));
        drawnow
    end   
end

if strcmp(Module,'OspreyLoadWater')
     if kk == 1    
        msg = sprintf('Loading raw data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
     msg = sprintf('Loading raw data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Loading raw data from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
%% Process
if strcmp(Module,'OspreyProcess') && ~MRSI
     if kk == 1    
        msg = sprintf('Processing data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
        msg = sprintf('Processing data from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Processing data from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
if strcmp(Module,'OspreyProcess') && MRSI && size(kk,2) > 1
     if kk(1) == 1   &&  kk(2) == 1
        msg = sprintf('Processing voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        fprintf(msg);
     else
     msg = sprintf('Processing voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Processing voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1)));
        drawnow
    end   
end
%% Fit
if strcmp(Module,'OspreyFit') && ~MRSI 
        msg = sprintf('\nFitting metabolite spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
    if GUI        
        set(progressText,'String' ,sprintf('Fitting metabolite spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
if strcmp(Module,'OspreyFit') && MRSI && size(kk,2) > 1
     if kk(1) == 1   &&  kk(2) == 1
        msg = sprintf('\nFitting metabolite spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        fprintf(msg);
     else
     msg = sprintf('\nFitting metabolite spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('\nFitting metabolite spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1)));
        drawnow
    end   
end

if strcmp(Module,'OspreyFitRef') && ~MRSI
        msg = sprintf('\nFitting water reference spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
    if GUI        
        set(progressText,'String' ,sprintf('Fitting water reference spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
if strcmp(Module,'OspreyFitRef') && MRSI && size(kk,2) > 1
     if kk(1) == 1   &&  kk(2) == 1
        msg = sprintf('\nFitting water reference spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        fprintf(msg);
     else
     msg = sprintf('\nFitting water reference spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('\nFitting water reference spectra from voxel %3i out of %3i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1)));
        drawnow
    end   
end

if strcmp(Module,'OspreyFitShortTE') && ~MRSI  
        msg = sprintf('\nFitting short-TE water spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
    if GUI        
        set(progressText,'String' ,sprintf('Fitting short-TE water spectra from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
if strcmp(Module,'OspreyFitShortTE') && MRSI && size(kk,2) > 1
     if kk(1) == 1   &&  kk(2) == 1
        msg = sprintf('\nFitting short-TE water spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        fprintf(msg);
     else
     msg = sprintf('\nFitting short-TE water spectra from voxel %5i out of %5i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1));
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('\nFitting short-TE water spectra from voxel %3i out of %3i from dataset %3i out of %3i total datasets...\n',kk(2),nDatasets(2), kk(1), nDatasets(1)));
        drawnow
    end   
end
%% Coreg
if strcmp(Module,'OspreyCoreg')
     if kk == 1    
        msg = sprintf('Coregistering voxel from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
     msg = sprintf('Coregistering voxel from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Coregistering voxel from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
%% Segmentation
if strcmp(Module,'OspreySeg')
     if kk == 1    
        msg = sprintf('\nSegmenting structural image from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
     msg = sprintf('\nSegmenting structural image from dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Coregistering voxel from dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
%% Quantification
if strcmp(Module,'OspreyQuant')
     if kk == 1    
        msg = sprintf('Quantifying dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        fprintf(msg);
     else
     msg = sprintf('Quantifying dataset %3i out of %3i total datasets...\n', kk, nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
     end
    if GUI        
        set(progressText,'String' ,sprintf('Quantifying dataset %3i out of %3i total datasets...\n', kk, nDatasets));
        drawnow
    end   
end
%% Done
if strcmp(Module,'done') && ~MRSI
    if GUI        
        set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',kk));
        pause(1);
    end
    fprintf('\n... done.\n Elapsed time %f seconds\n',kk);
end

if strcmp(Module,'MRSIdone')
    if GUI        
        set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',kk));
        pause(1);
    end
    fprintf('\n... done.\n Elapsed time %f seconds\n',kk);
end
if strcmp(Module,'Fulldone')
    if GUI        
        set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',kk));
        pause(1);
    end
    fprintf('\n... done.\n Elapsed time %f seconds\n',kk);
end
end