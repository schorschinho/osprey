function [MRSCont] = RunOspreyJob(jobFilePath)
%% [MRSCont] = RunOspreyJob(jobFilePath)
%   This function creates performs a full Osprey analysis on the supplied 
%   job file. This is the foundation for the compiled version of Osprey.
%
%   USAGE:
%       MRSCont = RunOspreyJob(jobFilePath)
%
%   INPUTS:
%       jobFilePath     = path to the Osprey jobfile (.m, .csv or .json).
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2022-05-17)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2022-05-17: First version of the code.

    MRSCont = OspreyJob(jobFilePath,0,'11');

    MRSCont = OspreyLoad(MRSCont);

    MRSCont = OspreyProcess(MRSCont);

    MRSCont = OspreyFit(MRSCont);

    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreyCoreg(MRSCont);
        if ~isempty(MRSCont.files_seg)
            seg = [num2str(~isempty(MRSCont.files_seg)) num2str(exist(MRSCont.files_seg{1}{1},'file'))];
        else
            seg = '00';
        end
        switch seg
            case '00' %Do normal segmentation (No files_seg)
                MRSCont = OspreySeg(MRSCont);
            case '11' %files_seg are defined and finished.
                MRSCont = OspreySeg(MRSCont);
            otherwise % (10) files_seg are defined, but not finished yet. (Do nothing)
        end
        
    end

    MRSCont = OspreyQuantify(MRSCont);

    MRSCont = OspreyOverview(MRSCont);

    for kk = 1 : MRSCont.nDatasets
        [MRSCont] = OspreyHTMLReport(MRSCont,kk);
    end

end

