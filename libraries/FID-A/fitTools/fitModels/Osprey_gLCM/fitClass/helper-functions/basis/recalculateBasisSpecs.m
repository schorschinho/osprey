function basisSet = recalculateBasisSpecs(basisSet)
%%  basisSet = recalculateBasisSpecs(basisSet)
% This function recalculates the basis spectra and ppm-axis of the basis set
%
%   USAGE:
%       basisSet = recalculateBasisSpecs(basisSet)
%
%   INPUTS:
%       basisSet        = Osprey basis set struct
%
%   OUTPUTS:
%       basisSet        = Osprey basis set struct with frequency domain      
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)   
%% Re-calcualte spectra and ppm axis

    basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);                   % Re-calculate spectra

    % Re-calculate ppm-axis
    f = [(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
    basisSet.ppm = f/(basisSet.Bo*42.577);                                  % Convert Hz to ppm
    basisSet.ppm = basisSet.ppm + basisSet.centerFreq;                      % Shift by center frequency
    
end