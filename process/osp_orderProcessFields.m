function [MRSCont] = osp_orderProcessFields(MRSCont)
%% [MRSCont] = osp_orderProcessFields(MRSCont)
%   This function orders the processed struct fields to match the GUI order
%
%   USAGE:
%       MRSCont = OspreyProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-08-19)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-08-19: First version of the code.

% Check that OspreyProcess has been run before
outputFolder = MRSCont.outputFolder;
fileID = fopen(fullfile(outputFolder, 'LogFile.txt'),'a+');
if ~MRSCont.flags.didProcess
    msg = 'Trying to sort processed struct fields, but data has not been processed yet. Run OspreyProcess first.';
    fprintf(fileID,msg);
    error(msg);
end

%% Define order based on the supplied files and the sequence type

if MRSCont.flags.isUnEdited || MRSCont.flags.isDWMRS
    numberstring = ['1' num2str(MRSCont.flags.hasMM) num2str(MRSCont.flags.hasRef) num2str(MRSCont.flags.hasWater)];
    switch numberstring
        case '1000'
            sortstring = {'A'}; 
        case '1100'
            sortstring = {'A','mm'};
        case '1001'
            sortstring = {'A','w'};
        case '1010'
            sortstring = {'A','ref'};
        case '1110'
            sortstring = {'A','mm','ref'};
        case '1101'
            sortstring = {'A','mm','w'};
        case '1111'
            sortstring = {'A','mm','ref','w'};       
        otherwise
            msg = 'Something is wrong in the processing!';
            fprintf(fileID,msg);
            error(msg);
    end
        
elseif MRSCont.flags.isMEGA
    numberstring = ['1' num2str(MRSCont.flags.hasRef) num2str(MRSCont.flags.hasWater)];
    switch numberstring
        case '100'
            sortstring = {'A','B','diff1','sum'}; 
        case '110'
            sortstring = {'A','B','diff1','sum','ref'};
        case '101'
            sortstring = {'A','B','diff1','sum','w'};
        case '111'
            sortstring = {'A','B','diff1','sum','ref','w'};
        otherwise
            msg = 'Something is wrong in the processing!';
            fprintf(fileID,msg);
            error(msg);
    end
    
elseif (MRSCont.flags.isHERMES) || (MRSCont.flags.isHERCULES)
     if ~(length(MRSCont.opts.editTarget) > 2)
        numberstring = ['1' num2str(MRSCont.flags.hasRef) num2str(MRSCont.flags.hasWater)];
        switch numberstring
            case '100'
                sortstring = {'A','B','C','D','diff1','diff2','sum'}; 
            case '110'
                sortstring = {'A','B','C','D','diff1','diff2','sum','ref'};
            case '101'
                sortstring = {'A','B','C','D','diff1','diff2','sum','w'};
            case '111'
                sortstring = {'A','B','C','D','diff1','diff2','sum','ref','w'};
            otherwise
                msg = 'Something is wrong in the processing!';
                fprintf(fileID,msg);
                error(msg);
        end  
     else
         numberstring = ['1' num2str(MRSCont.flags.hasRef) num2str(MRSCont.flags.hasWater)];
        switch numberstring
            case '100'
                sortstring = {'A','B','C','D','diff1','diff2','diff3','sum'}; 
            case '110'
                sortstring = {'A','B','C','D','diff1','diff2','diff3','sum','ref'};
            case '101'
                sortstring = {'A','B','C','D','diff1','diff2','diff3','sum','w'};
            case '111'
                sortstring = {'A','B','C','D','diff1','diff2','diff3','sum','ref','w'};
            otherwise
                msg = 'Something is wrong in the processing!';
                fprintf(fileID,msg);
                error(msg);
        end  
     end
else
    msg = 'No flag set for sequence type!';
    fprintf(fileID,msg);
    error(msg);
end

%% Sort the struct accordingly
if exist('sortstring','var') && (length(fieldnames(MRSCont.processed)) == length(sortstring))
    [MRSCont.processed, ~] = orderfields(MRSCont.processed,sortstring);
else
    msg = 'Something is wrong in the processing!';
    fprintf(fileID,msg);
    error(msg);    
end
fclose(fileID);
end