function [MRSCont] = osp_getRawBlocks(MRSCont,block_size)
%% [MRSCont] = osp_processMEGA(MRSCont, block_size)
%   This function performs the following steps to process MEGA-edited
%   (2-step) MRS data (e.g. MEGA-PRESS, MEGA-sLASER):
%       - Alignment of individual averages using robust spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Automated zero-order phase correction
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = osp_processMEGA(MRSCont, 'target');
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       block_size      = String. Can be 'GABA' or 'GSH'. Default: 'GABA'
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-22)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-05-27: First public version of the code.


MRSCont_bck = MRSCont;

if mod((MRSCont.raw{1}.rawAverages/2),block_size) == 0
    no_blocks = (MRSCont.raw{1}.rawAverages/2)/block_size;
    
    % Fake the blocks in as separate subjects for now... will do this in
    % additional dimensions later
    MRSCont.nDatasets = no_blocks;
    
    for nb = 1 : no_blocks
        MRSCont.raw{nb} =MRSCont_bck.raw{1};
        MRSCont.raw{nb}.rawAverages = block_size*2;
        MRSCont.raw{nb}.averages = block_size;
        MRSCont.raw{nb}.sz = [MRSCont_bck.raw{1}.sz(1) block_size MRSCont_bck.raw{1}.sz(3)];
        
        MRSCont.raw{nb}.fids = MRSCont_bck.raw{1}.fids(:,1+block_size*(nb-1):block_size*nb,:);
        MRSCont.raw{nb}.specs = MRSCont_bck.raw{1}.specs(:,1+block_size*(nb-1):block_size*nb,:);        
        MRSCont.files{nb} = MRSCont_bck.files{1};
        if MRSCont_bck.flags.hasRef
            MRSCont.files_ref{nb} = MRSCont_bck.files_ref{1};
            MRSCont.raw_ref{nb} = MRSCont_bck.raw_ref{1};
        end
        if MRSCont_bck.flags.hasWater
            MRSCont.files_w{nb} = MRSCont_bck.files_w{1};
            MRSCont.raw_w{nb} = MRSCont_bck.raw_w{1};
        end
        if ~isempty(MRSCont_bck.files_nii)
            MRSCont.files_nii{nb} = MRSCont_bck.files_nii{1};
        end
    end
    
else
    msg = ' rawAverages mod block_size is not zero.';
    fprintf(msg);
    error(msg);
end
end