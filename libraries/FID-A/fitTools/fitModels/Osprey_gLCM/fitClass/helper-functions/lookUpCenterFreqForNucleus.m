function centerFreq = lookUpCenterFreqForNucleus(nucleus)
% Returns the ppm value at the center of the spectrum depending
% on the nucleus.
if iscell(nucleus)
    nucleus = nucleus{1};
end
switch strtrim(nucleus)
    case '1H'
        centerFreq = 4.68;
    case '2H'
        centerFreq = 4.68;
    case '31P'
        centerFreq = 0;
    otherwise
        error('Nucleus %s not supported yet.', nucleus);
end

end