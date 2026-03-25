function [BASIS] = create_maximumEchoBasis(BASIS,tstart,side)

% Calculate time to add
tadd = BASIS.te-tstart;

tadd = tadd/1000;

pointsToAdd = tadd/BASIS.dwelltime;

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