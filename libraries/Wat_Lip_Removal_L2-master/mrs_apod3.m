function fids_apod = mrs_apod3( fids, BW, LB, Pos)
% MRS_APOD multiplies each FID by an exponential filter (a line broadening filter)  
% to eliminate sinc ringing caused by a truncated FID. 
% Exponential Filter: E(t)=exp(-pi*LB*t).  
% fids_apod = mrs_apod( fids, BW, LB )
% ARGS :
% fids = FIDs before apodization       (dim=[samples,avgs,dyns])
% BW = spectral bandwidth   (Hz)
% LB = the FWHM of the exponential filter (Hz)
% Pos = the spin echo position
% RETURNS:
% fids_apod = FIDs after apodization   (dim=[samples,avgs,dyns])
% EXAMPLE: 
% >> FID_apod = mrs_apod(FID, 4000, 5, 121); 
% >> figure; plot(real(FID));
% >> hold on; plot(real(FID_apod),'r');
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
% Copyright (c) 2013, University of Nottingham. All rights reserved.
% Edited by : Liangjie Lin
% PLACE  : Johns Hopkins University, Department of Radiology 

    dimnum = size(fids);
    xx =length(dimnum);
    % apodization
    if LB~=0 
        t = 0:(1/BW):((dimnum(3)-1)/BW);
        filter=exp(-pi*LB.*abs(t-((Pos-1)/BW)));        
        filter=repmat(filter,dimnum(2),1);
        
        if xx == 3
            for d=1:dimnum(1)
                fids_apod(d,:,:)=filter.*squeeze(fids(d,:,:));
            end
        end
        if xx == 4
            for d=1:dimnum(1)
                for e = 1:dimnum(4)
                    fids_apod(d,:,:,e)=filter.*squeeze(fids(d,:,:,e));
                end
            end
        end
            
    else
        fids_apod=fids;
    end     
end

