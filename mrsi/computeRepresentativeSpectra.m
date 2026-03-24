function MRSCont = computeRepresentativeSpectra(MRSCont,MetabSpecName)
%% MRSCont = computeRepresentativeSpectra(MRSCont)
%   Computes representative WM and GM
%   spectra according to Goryawala et al. MRM 2017 
%
%   USAGE:
%       MRSCont = computeRepresentativeSpectra(MRSCont)
%
%   INPUTS:
%       MRSCont    = Osprey data container.
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-10-31)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-10-31: First version of the code.



    %   D - Matrix of voxel spectra (n x m)
    %       n = number of voxels
    %       m = number of spectral points
    %   W - Fractional WM and GM content (n x 2)

    

    % Check that OspreyCoreg has been run before
    if ~MRSCont.flags.didSeg
        error('Trying to calculate representatitve spectra, but segmentation has not been performed yet. Run OspreySeg first.')
    end


    for kk = 1:MRSCont.nDatasets
        % Lets get the data we need by first reshaping it 
        D = reshape(MRSCont.processed.(MetabSpecName){kk}.specs(:,:,:,MRSCont.opts.MRSI.RepSpectra.SliceIndices),MRSCont.processed.(MetabSpecName){kk}.sz(1),[]);

        WM = reshape(squeeze(MRSCont.seg.tissue.fWM(kk,:,:,MRSCont.opts.MRSI.RepSpectra.SliceIndices)),1,[]);
        GM = reshape(squeeze(MRSCont.seg.tissue.fGM(kk,:,:,MRSCont.opts.MRSI.RepSpectra.SliceIndices)),1,[]); 
        tWM = WM;
        tGM = GM;

        D = D(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)';
        WM = WM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)'./(tWM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)' + tGM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)');
        GM = GM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)'./(tWM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)' + tGM(:,(tWM + tGM) > MRSCont.opts.MRSI.RepSpectra.fGMpfWM)');
             
        W = cat(2,WM,GM);
        
        % Compute representative spectra using least squares solution
        S = pinv(W) * D;
        S = S';



        MRSCont.RepSpectra = MRSCont.processed.(MetabSpecName);
        MRSCont.RepSpectra{kk}.specs = S;
        MRSCont.RepSpectra{kk}.fids=ifft(fftshift(S,MRSCont.RepSpectra{kk}.dims.t),[],MRSCont.RepSpectra{kk}.dims.t);
        MRSCont.RepSpectra{kk}.sz = size(S);
        MRSCont.RepSpectra{kk}.dims.extras = 1;
        MRSCont.RepSpectra{kk}.dims.Xvoxels = 0;
        MRSCont.RepSpectra{kk}.dims.Yvoxels = 0;
        MRSCont.RepSpectra{kk}.dims.Zvoxels = 0;

    end
    
end