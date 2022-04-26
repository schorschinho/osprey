% fit_makeMMBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [BASIS] = fit_makeBasis(folder, addMMFlag, sequence, editTarget)
% 
% DESCRIPTION:
% Generates a basis set in FID-A structure. The code will search all *.mat
% files in the input folder for FID-A structures with simulated spectra. It
% also performs sanity checks on the simulation parameters, and returns
% warnings if parameters are not identical for all parameters.
% 
% INPUTS:
% folder    = folder containing *.mat files representing FID-A structures
% addMMFlag = Flag to decide whether MM and lipid basis functions should be
%               added to the basis set.
%             OPTIONS:  1 = Add MM+lip (Default)
%                       0 = Don't add MM+lip
% sequence  = sequence type
%             OPTIONS:  'unedited' (default)
%                       'MEGA'
%                       'HERMES'
%                       'HERCULES'
% editTarget= Target molecule of edited data.
%             OPTIONS:  'GABA'
%                       'GSH'
%                       '

%
% OUTPUTS:
% BASIS     = Simulated basis set in FID-A structure format. 

function [BASIS] = fit_makeMMBasis(inBASIS)

    BASIS = inBASIS;
    
    BASIS.fids  = BASIS.fids .* BASIS.scale;
    BASIS.specs = BASIS.specs .* BASIS.scale;
    
    
    n = BASIS.n;
    sw = BASIS.spectralwidth;
    Bo = BASIS.Bo;
    centerFreq = BASIS.centerFreq;
    % The amplitude and FWHM values are determined as for the LCModel and
    % TARQUIN algorithms (see Wilson et al., MRM 2011).
    hzppm = Bo*42.577;
    
    % To scale the amplitudes correctly, we first need to determine the
    % area of the 3.027 ppm CH3 signal of creatine
    [CrArea] = detCrArea(BASIS);
    oneProtonArea = CrArea/3;
    
    % Next, we determine the area of a Gaussian singlet with nominal area 1
    testGaussian    = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,0,1);
    testGaussian    = op_dccorr(testGaussian,'p');
    gaussianArea    = sum(real(testGaussian.specs));
    
    % Now we know the scaling factor to generate MM/lipid signals with the
    % correct relative scaling with respect to the CH3 signal
    MM09            = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,0.91,3*oneProtonArea/gaussianArea);
    MMBase.MM09     = op_dccorr(MM09,'p');
    MM12            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.21,2*oneProtonArea/gaussianArea);
    MMBase.MM12     = op_dccorr(MM12,'p');
    MM15            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.51,2*oneProtonArea/gaussianArea);
    MMBase.MM15     = op_dccorr(MM15,'p');
    MMBase.MM15.fids = -MMBase.MM15.fids;
    MMBase.MM15.specs = -MMBase.MM15.specs;
    MM18            = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,1.8,2*oneProtonArea/gaussianArea);
    MMBase.MM18     = op_dccorr(MM18,'p');
    MM20           = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,2.08,3*oneProtonArea/gaussianArea);
    MMBase.MM20     = op_dccorr(MM20,'p');
    MMBase.MM20.fids = -MMBase.MM20.fids;
    MMBase.MM20.specs = -MMBase.MM20.specs;
    MM23           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.3,0.33*oneProtonArea/gaussianArea);
    MMBase.MM23     = op_dccorr(MM23,'p');
    MMBase.MM23.fids = MMBase.MM23.fids;
    MMBase.MM23.specs = MMBase.MM23.specs;
    MM3co  = op_gaussianPeak(n,sw,Bo,centerFreq,14,3,2*oneProtonArea/gaussianArea);
    MM3co  = op_dccorr(MM3co,'p');
    MMBase.MM3co     = op_dccorr(MM3co,'p');
    
%     MM20a           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.08,1.33*oneProtonArea/gaussianArea);
%     MM20b           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.25,0.33*oneProtonArea/gaussianArea);
%     MM20c           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.95,0.33*oneProtonArea/gaussianArea);
%     MM20d           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,3.0,0.4*oneProtonArea/gaussianArea);
%     MM20            = op_addScans(MM20a,MM20b); MM20 = op_addScans(MM20,MM20c); MM20 = op_addScans(MM20,MM20d);
%     MMBase.MM20     = op_dccorr(MM20,'p');

    MMLips = {'MM09','MM12','MM15','MM18','MM20','MM23','MM3co'};

    % Now copy over the names, fids, and specs into the basis set structure
    BASIS = rmfield(BASIS,'fids');
    BASIS = rmfield(BASIS,'specs');
    BASIS = rmfield(BASIS,'name');
    BASIS.fids(:,1)   = MMBase.(MMLips{1}).fids;
    BASIS.specs(:,1)  = MMBase.(MMLips{1}).specs;
    BASIS.name{1}       = MMLips{1};
    for rr = 2:length(MMLips)
        BASIS.name{rr}       = MMLips{rr};
        BASIS.fids(:,rr)   = MMBase.(MMLips{rr}).fids;
        BASIS.specs(:,rr)  = MMBase.(MMLips{rr}).specs;
    end
 

% Copy over the FID, specs, dims, and the metabolite names
BASIS.flags.addedMM     = 1;
BASIS.nMM               = length(MMLips);
BASIS.nMets            = 0;
BASIS.sz                = size(BASIS.fids);

BASIS.fids  = BASIS.fids ./ BASIS.scale;
BASIS.specs = BASIS.specs ./ BASIS.scale;


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
refWindow = in.specs(ppm>ppmmin & ppm<ppmmax, idx_Cr,1);
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