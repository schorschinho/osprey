function [BASISdiff] = osp_addDiffMMPeaks(BASISdiff,BASISoff,fitOpts)
    n = BASISdiff.n;
    sw = BASISdiff.spectralwidth;
    Bo = BASISdiff.Bo;
    centerFreq = BASISdiff.centerFreq;
    % The amplitude and FWHM values are determined as for the LCModel and
    % TARQUIN algorithms (see Wilson et al., MRM 2011).
    hzppm = Bo*42.577;

    % To scale the amplitudes correctly, we first need to determine the
    % area of the 3.027 ppm CH3 signal of creatine
    [CrArea] = detCrArea(BASISoff);
    oneProtonArea = CrArea/3;
    
    % Next, we determine the area of a Gaussian singlet with nominal area 1
    testGaussian    = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,0,1);
    testGaussian    = op_dccorr(testGaussian,'p');
    gaussianArea    = sum(real(testGaussian.specs));
    
     % Find the GABA basis function
    idx_GABA          = find(strcmp(BASISdiff.name,'GABA'));
    if isempty(idx_GABA)
        error('No basis function with nametag ''GABA'' found! Abort!');
    end
 
    if strcmp(fitOpts.coMM3, '1to1GABA') % 1:1 GABA and co-edited MM3 model
        MM3co  = op_gaussianPeak(n,sw,Bo,centerFreq,fitOpts.FWHMcoMM3,3,2*oneProtonArea/gaussianArea);
        MM3co  = op_dccorr(MM3co,'p');
        BASISdiff.fids(:,idx_GABA)  = BASISdiff.fids(:,idx_GABA) +  MM3co.fids;
        BASISdiff.specs(:,idx_GABA) = BASISdiff.specs(:,idx_GABA) +  MM3co.specs;
    end    
    
    if strcmp(fitOpts.coMM3, '3to2MM') % 3:2 MM09 and co-edited MM3 model
        idx_MM          = find(strcmp(BASISdiff.name,'MM09'));
        if isempty(idx_MM)
            error('No basis function with nametag ''MM09'' found! Abort!');
        end
        MM3co  = op_gaussianPeak(n,sw,Bo,centerFreq,fitOpts.FWHMcoMM3,3,2*oneProtonArea/gaussianArea);
        MM3co  = op_dccorr(MM3co,'p');
        BASISdiff.fids(:,idx_MM)  = BASISdiff.fids(:,idx_MM) +  MM3co.fids;
        BASISdiff.specs(:,idx_MM) = BASISdiff.specs(:,idx_MM) +  MM3co.specs;
    end   
    
    if strcmp(fitOpts.coMM3, '3to2MMsoft') % 3:2 MM09 and co-edited MM3 model
        idx_MM          = find(strcmp(BASISdiff.name,'MM09'));
        if isempty(idx_MM)
            error('No basis function with nametag ''MM09'' found! Abort!');
        end
        MM3co  = op_gaussianPeak(n,sw,Bo,centerFreq,fitOpts.FWHMcoMM3,3,2*oneProtonArea/gaussianArea);
        MM3co  = op_dccorr(MM3co,'p');
        BASISdiff.name{BASISdiff.nMets+BASISdiff.nMM+1}       = 'MM3co';
        BASISdiff.fids(:,BASISdiff.nMets+BASISdiff.nMM+1)   = MM3co.fids;
        BASISdiff.specs(:,BASISdiff.nMets+BASISdiff.nMM+1)  = MM3co.specs;    
        BASISdiff.nMM = BASISdiff.nMM + 1;
    end 

    if strcmp(fitOpts.coMM3, 'freeGauss') || strcmp(fitOpts.coMM3, 'fixedGauss') || strcmp(fitOpts.coMM3, '1to1GABAsoft')...
        || strcmp(fitOpts.coMM3, '2to3GABAsoft') || strcmp(fitOpts.coMM3, '3to2GABAsoft')% free co-edited MM3 model and soft constraints model
        MM3co  = op_gaussianPeak(n,sw,Bo,centerFreq,fitOpts.FWHMcoMM3,3,2*oneProtonArea/gaussianArea);
        MM3co  = op_dccorr(MM3co,'p');
        BASISdiff.name{BASISdiff.nMets+BASISdiff.nMM+1}       = 'MM3co';
        BASISdiff.fids(:,BASISdiff.nMets+BASISdiff.nMM+1)   = MM3co.fids;
        BASISdiff.specs(:,BASISdiff.nMets+BASISdiff.nMM+1)  = MM3co.specs;    
        BASISdiff.nMM = BASISdiff.nMM + 1;
    end
    
    if strcmp(fitOpts.coMM3, 'freeGauss')
        BASISdiff.coMM3.MMArea = 2*oneProtonArea/gaussianArea;
        BASISdiff.coMM3.n = n;
        BASISdiff.coMM3.sw = sw;
        BASISdiff.coMM3.Bo = Bo;
        BASISdiff.coMM3.centerFreq = centerFreq;        
    end
    
end

% detCrArea.m
% Georg Oeltzschner, Johns Hopkins University 2020
% 
% USAGE:
% [CrArea] = detCrArea(buffer);
% 
% DESCRIPTION:
% Finds the creatine spectrum in the temporary basis set buffer, then fits
% a Lorentzian to the 3.027 ppm CH3 creatine singlet to determine its area.
% Subsequently, macromolecule and lipid basis functions are scaled
% accordingly.
% 
% INPUTS:
% in        = a temporary buffer containing simulated basis functions
%
% OUTPUTS:
% CrArea    = Estimated area under the 3.027 ppm CH3 Cr singlet.


function [CrArea] = detCrArea(in);

% Find the creatine basis function
idx_Cr          = find(strcmp(in.name,'Cr'));
if isempty(idx_Cr)
    error('No basis function with nametag ''Cr'' found! Abort!');
end

%[~, idx_3027]   = min(abs(buffer.ppm(:,1)-3.027));

% Determine the window where we are going to look for the peak.
ppm = in.ppm(1,:);
ppmmin = 3.027 - 0.4;
ppmmax = 3.027 + 0.4;
refWindow = in.specs(ppm>ppmmin & ppm<ppmmax, idx_Cr);
ppmWindow = in.ppm(ppm>ppmmin & ppm<ppmmax);

% Find the maximum and its index
maxRef_index    = find(abs(real(refWindow)) == max(abs(real((refWindow)))));
maxRef          = real(refWindow(maxRef_index));

% Determine an initial estimate for the FWHM
% Peak lines can be super narrow, so overestimate it slightly
gtHalfMax   = find(abs(real(refWindow)) >= 0.4*abs(maxRef));
FWHM1       = abs(ppmWindow(gtHalfMax(1)) - ppmWindow(gtHalfMax(end)));
FWHM1       = FWHM1*(42.577*in.Bo(1));  %Assumes proton.

% Determine an initial estimate for the center frequency of the Cr peak
crFreq = ppmWindow(maxRef_index);

% Set up the fit
parsGuess=zeros(1,5);
parsGuess(1) = maxRef;  % amplitude
parsGuess(2) = (5*in.Bo/3)/(42.577*in.Bo); %FWHM.  Assumes Proton.  LW = 5/3 Hz/T.   % FWHM. Assumes Proton.
parsGuess(3) = crFreq;  % center frequency
parsGuess(4) = 0;       % baseline offset
parsGuess(5) = 0;       % phase
    
% Run first guess
yGuess  = op_lorentz(parsGuess, ppmWindow);
parsFit = nlinfit(ppmWindow, real(refWindow'), @op_lorentz, parsGuess);
yFit    = op_lorentz(parsFit, ppmWindow);
    
% figure;
% plot(ppmWindow,refWindow,'.',ppmWindow,yGuess,':',ppmWindow,yFit);
% legend('data','guess','fit');

CrArea = sum(yFit);

end