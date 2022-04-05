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

numberstring = ['1' num2str(MRSCont.flags.hasMM) num2str(MRSCont.flags.hasRef) num2str(MRSCont.flags.hasWater) num2str(MRSCont.flags.hasMMRef)];
switch numberstring
    case '10000'
        sortstring = {'metab'};
    case '10110'
        sortstring = {'metab','ref','w'};      
    case '10100'
        sortstring = {'metab','ref'};  
    case '10010'
        sortstring = {'metab','w'};          
    case '11000'
        sortstring = {'metab','mm'}; 
    case '11001'
        sortstring = {'metab','mm','mm_ref'};         
    case '11101'
        sortstring = {'metab','mm','ref','mm_ref'};         
    case '11110'
        sortstring = {'metab','mm','ref','w'}; 
    case '11011'
        sortstring = {'metab','mm','w','mm_ref'};          
    case '11111'
        sortstring = {'metab','mm','ref','w','mm_ref'};       
    otherwise
        msg = 'Something is wrong in the processing!';
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