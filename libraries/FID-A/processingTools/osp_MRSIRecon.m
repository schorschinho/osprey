function [MRSCont] = osp_MRSIRecon(MRSCont)
%% [MRSCont] = osp_MRSIRecon(MRSCont)
%   This function reads MRSI data.
%
%   USAGE:
%       [MRSCont] = osp_MRSIRecon(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2021-01-04)
%       hzoelln2@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-01-04: First version of the code.

% Close any remaining open figures
close all;
if (MRSCont.flags.isMRSI==0)
    error('ERROR:  This is not a MRSI dataset!  Aborting!');
end
% fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');


%% Get the data (loop over all datasets)
% refLoadTime = tic;
% reverseStr = '';
% if MRSCont.flags.isGUI
%     progressText = MRSCont.flags.inProgress;
% end

for kk = 1:MRSCont.nDatasets

%     msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
%     fprintf([reverseStr, msg]);
%     fprintf(fileID,[reverseStr, msg]);
%     reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     if MRSCont.flags.isGUI        
%         set(progressText,'String' ,sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
%     end
            
    switch MRSCont.vendor
        case 'Philips'
            switch MRSCont.datatype
                case 'SDAT'
                    raw = op_sortMRSIsdat(MRSCont.raw{kk});
                    MRSCont.raw{kk}      = raw;
                    if MRSCont.flags.hasWater
                        raw_w = op_sortMRSIsdat(MRSCont.raw_w{kk});
                        MRSCont.raw_w{kk}      = raw_w;
                    end
                otherwise
                    msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                    fprintf(fileID,msg);
                    error(msg);
            end
        otherwise
            msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
            fprintf(fileID,msg);
            error(msg);
    end
end
% fprintf('... done.\n');
% time = toc(refLoadTime);
% if MRSCont.flags.isGUI        
%     set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
%     pause(1);
% end
% fprintf(fileID,'... done.\n Elapsed time %f seconds\n',time);
% fclose(fileID);
% MRSCont.runtime.Load = time;
end