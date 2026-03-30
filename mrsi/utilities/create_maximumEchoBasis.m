function [BASIS] = create_maximumEchoBasis(BASIS,tstart,side)
%% [BASIS] = create_maximumEchoBasis(BASIS,tstart,side)
%   Creates a maximum echo basisset from a normal Osprey basisset
%
%   USAGE:
%       [BASIS] = create_maximumEchoBasis(BASIS,tstart,side)
%
%   INPUTS:
%       BASIS  = Osprey basis set.
%       tstart = When did the ADC start (in ms)
%       side   = What do you want to export see options below
%
%   OUTPUTS:
%       BASIS  = Osprey basis set.
%
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.
%% Create basis
% Calculate time to add
tadd = BASIS.te-tstart;
tadd = tadd/1000;
pointsToAdd = tadd/BASIS.dwelltime;

% Different export options
switch side
    case 'max'
        for ii = 1 : BASIS.nMets
            fids(:,ii) = [conj(flipud(BASIS.fids(2:pointsToAdd,ii)));BASIS.fids(1:end-(pointsToAdd-1),ii)];
        end  
    case 'leftflipconj' 
        fids = zeros(size(BASIS.fids));
        for ii = 1 : BASIS.nMets
            fids(1:pointsToAdd,ii,:) = conj(flipud(BASIS.fids(1:pointsToAdd,ii,:)));
        end 
    case 'leftflip' 
        fids = zeros(size(BASIS.fids));
        for ii = 1 : BASIS.nMets
            fids(1:pointsToAdd,ii) = ((BASIS.fids(1:pointsToAdd,ii)));
        end 
end

specs = fftshift(fft(fids, [], 1), 1);
BASIS.fids = fids;
BASIS.specs = specs;
end