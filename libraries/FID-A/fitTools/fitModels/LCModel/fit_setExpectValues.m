% fit_setExpectValues.m
% Georg Oeltzschner, Johns Hopkins University 2020.
%
% USAGE:
% [EXT2, SDT2, SDSH] = fit_setExpectValues(dataToFit, basisSet)
% 
% DESCRIPTION:
% This function determines the expectation and standard deviation values 
% for the delta (1/T2) parameters describing the Lorentzian linebroadening 
% applied to each basis function in the LCModel algorithm.
%
% These values are used in a penalty term to regularize the linebroadening
% that is applied across the basis functions, and the frequency shifts
% that they can have. (Eq. (4) in the original LCModel paper, 
% Provencher, Magn Reson Med 30:672-679 (1993)).
%
% 
% OUTPUTS:
% EXT2        = Vector of delta (1/T2) expectation values [Hz]
% SDT2        = Vector of delta (1/T2) standard deviation values [Hz]
% SDSH        = Vector of frequency shift standard deviation values [Hz]
%
% INPUTS:
% dataToFit   = FID-A data structure with the dataset to be fit
% basisSet    = FID-A basis set container

function [EXT2, SDT2, SDSH] = fit_setExpectValues(dataToFit, basisSet)

%%% 1. GET DEFAULT VALUES %%%
% See LCModel manual, Section 11.14 (p. 149)

% First, a scaling factor to roughly account for the field strength
% dependence of T2 is determined
hzpppm = dataToFit.txfrq*1e-6;
scalingT2 = sqrt(hzpppm / 85.15); % scaling factor to account for T2 decrease with field strength

% There are some default values listed for the parameters
defEXT2 = 2.0;     % [Hz] - expectation value for Lorentzian linebroadening
defSDT2 = 0.4;     % [Hz] - standard deviation for Lorentzian linebroadening
defSDSH = 0.004*dataToFit.txfrq*1e-6;   % [Hz] - standard deviation for frequency shifts


%%% 2. LOOP OVER BASIS FUNCTIONS AND MAKE ADJUSTMENTS %%%
% NAA and NAAG, by default, get a lower SDSH to improve their separation
% Extract the names of the basis functions
basisNames = basisSet.name;

% Set up the initial output vectors
EXT2 = zeros(size(basisNames));
SDT2 = zeros(size(basisNames));
SDSH = zeros(size(basisNames));

% First, loop over the metabolite basis functions
nMets = basisSet.nMets;
nMM   = basisSet.nMM;
for rr = 1:nMets
    % First, do the EXT2
    EXT2(rr) = defEXT2 * scalingT2;
    SDT2(rr) = defSDT2 * scalingT2;
    
    % Then, do the SDSH
    % NAA and NAAG get less 'wiggle room' to improve their separation
    if strcmp(basisNames{rr}, 'NAA') || strcmp(basisNames{rr}, 'NAAG')
        SDSH(rr) = 0.002*hzpppm;
    else
        SDSH(rr) = defSDSH;
    end
end

% Next, loop over the MM basis functions.
% Provide the default expectation values for the frequency shift
% The EXT2 and SDT2 calculation is a bit more complicated - need to figure
% out exactly how it's done!
if nMM > 0
    for rr = 1:nMM
        if strcmp(basisNames{rr+nMets}, 'MM09')
            gaussFWHM(rr)   = 0.14 * hzpppm;
            targetFWHM(rr)  = 0.17 * hzpppm; 
            SDtargetFWHM(rr)  = 0.015 * hzpppm;
            SDSH(rr+nMets)        = 0.02 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'MM12')
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.20 * hzpppm;
            SDtargetFWHM(rr)  = 0.02 * hzpppm;
            SDSH(rr+nMets)        = 0.01 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'MM14')
            gaussFWHM(rr)   = 0.17 * hzpppm;
            targetFWHM(rr)  = 0.20 * hzpppm;
            SDtargetFWHM(rr)  = 0.02 * hzpppm;
            SDSH(rr+nMets)        = 0.02 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'MM17')
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.17 * hzpppm;
            SDtargetFWHM(rr)  = 0.02 * hzpppm;
            SDSH(rr+nMets)        = 0.03 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'MM20')
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.18 * hzpppm;
            SDtargetFWHM(rr)  = 0.01 * hzpppm;
            SDSH(rr+nMets)        = 0.005* hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'Lip09')
            gaussFWHM(rr)   = 0.14 * hzpppm;
            targetFWHM(rr)  = 0.19 * hzpppm;
            SDtargetFWHM(rr)  = 0.035 * hzpppm;
            SDSH(rr+nMets)        = 0.02 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'Lip13')
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.20 * hzpppm;
            SDtargetFWHM(rr)  = 0.035 * hzpppm;
            SDSH(rr+nMets)        = 0.01 * hzpppm;
        elseif strcmp(basisNames{rr+nMets}, 'Lip20')
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.20 * hzpppm;
            SDtargetFWHM(rr)  = 0.025 * hzpppm;
            SDSH(rr+nMets)        = 0.01 * hzpppm;
        else
            gaussFWHM(rr)   = 0.15 * hzpppm;
            targetFWHM(rr)  = 0.20 * hzpppm;
            SDtargetFWHM(rr)  = 0.02 * hzpppm;
            SDSH(rr+nMets)        = 0.01 * hzpppm;
        end
        
        % For now, set them all to defEXT2;
        % First, do the EXT2
        EXT2(rr+nMets) = defEXT2 * scalingT2;
        SDT2(rr+nMets) = defSDT2 * scalingT2;
      
    end
    
%     % Now, try to derive proper EXT2 values
%     for rr = 1:length(gaussFWHM)
%         f_g = gaussFWHM(rr);
%         f_l = linspace(0, 20, 20000);
%         % The linewidth of a Voigt profile can be approximated by the
%         % following formula, found on:
%         % https://en.wikipedia.org/wiki/Voigt_profile
%         % (f_l = Lorentzian FWHM, f_g = Gaussian FWHM, f_v = Voigtian FWHM)
%         f_v = 0.5346*f_l + sqrt(0.2166*f_l.^2 + f_g^2);
% 
%         [~, idx_optimalLorentzFWHM] = min(abs(f_v - targetFWHM(rr)));
%         lorentzFWHM(rr) = f_l(idx_optimalLorentzFWHM);
%         % R is the 1/T2 parameter; so we get R from the fundamental
%         % Lorentzian relationship FWHM = R/pi;
%         EXT2(rr+nMets) = lorentzFWHM(rr) * pi;
%         EXT2(rr+nMets) = EXT2(rr+nMets) * scalingT2;
%         
%         % We obtain the standard deviations just by relating them linearly
%         % to the expectation values
%         EXtoSD(rr) = SDtargetFWHM(rr)/targetFWHM(rr);
%         SDT2(rr+nMets) = EXT2(rr+nMets) * EXtoSD(rr) * scalingT2;
%     end


end


end