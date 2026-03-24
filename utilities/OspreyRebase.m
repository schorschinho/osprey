function [MRSCont] = OspreyRebase(MRSCont)
%% [MRSCont] = OspreyRebase(MRSCont)
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

oldOutputFolder = MRSCont.outputFolder;

outputFolder = uipickfiles('FilterSpec',userpath,'REFilter', '\','NumFiles',1,'Prompt','Select the derivative folder in the new location)');
outputFolder = outputFolder{1};
MRSCont.outputFolder = outputFolder;
outputFile = MRSCont.outputFile;

jobFolder = uipickfiles('FilterSpec',MRSCont.outputFolder,'REFilter', '\','NumFiles',1,'Prompt','Select the folder containing the jobFile in the new location)');
jobFolder = jobFolder{1};

[oldPath,~,~] = fileparts(MRSCont.loadedJob);
MRSCont.loadedJob = strrep(MRSCont.loadedJob,oldPath,jobFolder);

dataFolder = uipickfiles('FilterSpec',MRSCont.outputFolder,'REFilter', '\','NumFiles',1,'Prompt','Select the data folder in the new location)');
dataFolder = dataFolder{1};

prompt = 'Input the path to the folder contianing all subject folders on the server: \n ';
oldPath = input(prompt)

if MRSCont.flags.hasFiles
    for kk = 1 : MRSCont.nDatasets
        MRSCont.files{kk} = strrep(MRSCont.files{kk},oldPath,dataFolder);
    end
end

if MRSCont.flags.hasRef
   for kk = 1 : MRSCont.nDatasets
        MRSCont.files_ref{kk} = strrep(MRSCont.files_ref{kk},oldPath,dataFolder);
    end 
end

if MRSCont.flags.hasWater
   for kk = 1 : MRSCont.nDatasets
        MRSCont.files_w{kk} = strrep(MRSCont.files_w{kk},oldPath,dataFolder);
    end 
end

if MRSCont.flags.hasMM
   for kk = 1 : MRSCont.nDatasets
        MRSCont.files_mm{kk} = strrep(MRSCont.files_mm{kk},oldPath,dataFolder);
    end 
end

if ~isempty(MRSCont.files_nii)
   for kk = 1 : MRSCont.nDatasets
        MRSCont.files_nii{kk} = strrep(MRSCont.files_nii{kk},oldPath,dataFolder);
    end 
end

if MRSCont.flags.didCoreg
   for kk = 1 : MRSCont.nDatasets
       MRSCont.coreg.vol_image{kk}.fname  = MRSCont.files_nii{kk};
       MRSCont.coreg.vol_mask{kk}.fname  = strrep(MRSCont.coreg.vol_mask{kk}.fname,oldOutputFolder,outputFolder);
   end
end

% Save everyting back
if  MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end