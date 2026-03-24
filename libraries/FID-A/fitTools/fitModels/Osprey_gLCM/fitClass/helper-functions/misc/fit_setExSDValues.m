% fit_setExSDValues.m
% Georg Oeltzschner, Johns Hopkins University 2020.
%
% USAGE:
% [EXT2, SDT2, SDSH] = fit_setExSDValues(dataToFit, basisSet, MMLipConfig, gauss_rt2)
% 
% DESCRIPTION:
% This function determines the expectation (EXT2) and standard deviation 
% (SDT2) values for the delta (1/T2) parameters describing the Lorentzian 
% linebroadening applied to each basis function in the LCModel algorithm, 
% as well as the standard deviation of the frequency shifts (SDSH).
%
% It uses an LCModel-type definition of MM and lipid signals to determine
% these values for internally simulated MMs and lipids.
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
% dataToModel = FID-A data structure with the dataset to be fit
% basisSet    = FID-A basis set container
% MMLipConfig = structure generated from an MM/Lipid definition JSON
%               file
% gauss_rt2   = Flag to determine the method for calculating the EXT2 and
%               SDT2 values (see LCModel v6.3 source code, line 4199f).
%               If set to 0 (default), these values will be determined
%               *correctly*; if set to 1, they will be determined *wrongly*.
%               This wrong way is, curiously, the LCModel default to 
%               maintain backwards compatibility with older LCModel 
%               versions, although it leads, as the source code states, to 
%               excessive linebroadening of the MM functions/
%               (Default: 0 or false)

function [EXT2, SDT2, SDSH] = fit_setExSDValues(dataToModel, basisSet, MMLipConfig, gauss_rt2)

%%% 0. PARSE INPUT %%%
if nargin < 4
    gauss_rt2 = 0;
    if nargin < 3
        MMLipConfig = false;
    end
end

%%% 1. GET DEFAULT VALUES %%%
% See LCModel manual, Section 11.14 (p. 149)

% First, a scaling factor to roughly account for the field strength
% dependence of T2 is determined
hzpppm = dataToModel.txfrq*1e-6;
scalingT2 = sqrt(hzpppm / 85.15); % scaling factor to account for T2 decrease with field strength

% There are some default values listed for the parameters
defEXT2 = 2.0;     % [Hz] - expectation value for Lorentzian linebroadening
defSDT2 = 0.4;     % [Hz] - standard deviation for Lorentzian linebroadening
defSDSH = 0.004*hzpppm;   % [Hz] - standard deviation for frequency shifts


%%% 2. LOOP OVER METABOLITE BASIS FUNCTIONS AND MAKE ADJUSTMENTS %%%
% NAA and NAAG, by default, get a lower SDSH to improve their separation
% Extract the names of the basis functions
basisNames = basisSet.name;

% Set up the initial output vectors
EXT2 = zeros(size(basisNames));
SDT2 = zeros(size(basisNames));
SDSH = zeros(size(basisNames));

% First, loop over the metabolite basis functions
nMets = basisSet.nMets;
for rr = 1:nMets
    % First, do the EXT2
    EXT2(rr) = defEXT2 * scalingT2;
    SDT2(rr) = defSDT2 * scalingT2;
    
    % Then, do the SDSH
    % NAA and NAAG get less 'wiggle room' to improve their separation
    if strcmp(basisNames{rr}, 'NAA') || strcmp(basisNames{rr}, 'NAAG')
        SDSH(rr) = 0.002*hzpppm;
    else if strcmp(basisNames{rr}, 'GPC') || strcmp(basisNames{rr}, 'PCh') || ...
            strcmp(basisNames{rr}, 'Cho') || strcmp(basisNames{rr}, 'fCho')
            SDSH(rr) = 0.006*hzpppm;
    else
            SDSH(rr) = defSDSH;
        end
    end
end

%%% 2. LOOP OVER MM BASIS FUNCTIONS AND MAKE ADJUSTMENTS %%%
% Next, loop over the MM basis functions.
if isstruct(MMLipConfig)
    nMM   = basisSet.nMM;
    if nMM > 0
        for rr = 1:nMM
            % Find the name of the MM basis function in the MM/lipid config
            % file and extract the parameterizations
            namesMMLipConfig = fieldnames(MMLipConfig);
            whichMMIdx = find(strcmp(fieldnames(MMLipConfig),basisNames{rr+nMets}));
            if whichMMIdx
                sifwmn(rr)      = MMLipConfig.(namesMMLipConfig{whichMMIdx}).fwmin(1) * hzpppm;
                sifwex(rr)      = MMLipConfig.(namesMMLipConfig{whichMMIdx}).fwex * hzpppm;
                sifwsd(rr)      = MMLipConfig.(namesMMLipConfig{whichMMIdx}).sdfw * hzpppm;
                SDSH(rr+nMets)  = MMLipConfig.(namesMMLipConfig{whichMMIdx}).sdppm * hzpppm;
            else
                % Some defaults
            end

        end

        % Determine the EXT2 and SDT2 values
        for rr = 1:nMM
            % Calculate EXT2
            % LCModel line 4205ff.
            if gauss_rt2
                fwhm_ex(rr) = sqrt(sifwex(rr)^2 - sifwmn(rr)^2);
            else
                fwhm_ex(rr) = sifwex(rr) - sifwmn(rr);
            end
            EXT2(rr+nMets) = fwhm_ex(rr) * pi;
            
            % Calculate SDT2
            % LCModel line 4217ff.
            if gauss_rt2
                SDT2(rr+nMets) = pi * sqrt((sifwex(rr) + sifwsd(rr))^2 - sifwmn(rr)^2) - EXT2(rr+nMets);
            else
                SDT2(rr+nMets) = pi * sifwsd(rr);
            end

        end

    end

end


end