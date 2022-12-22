function [basisSet]=osp_recalculate_basis_specs(basisSet)
    % This function recalculates the basis spectra and ppm-axis of the
    % basis set

    basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

    % Calcualte ppm-axis
    f = [(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
    basisSet.ppm = f/(basisSet.Bo*42.577);
    basisSet.ppm=basisSet.ppm + basisSet.centerFreq;
end