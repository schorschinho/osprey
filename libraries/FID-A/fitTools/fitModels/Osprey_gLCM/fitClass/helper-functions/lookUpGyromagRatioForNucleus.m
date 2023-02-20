function gamma = lookUpGyromagRatioForNucleus(nucleus)
% Returns the gyromagnetic ratio [MHz/T] depending
% on the nucleus.
if iscell(nucleus)
    nucleus = nucleus{1};
end
switch strtrim(nucleus)
    case '1H'
        gamma = 42.577478518;
    case '2H'
        gamma = 6.536;
    case '13C'
        gamma = 10.7084;
    case '19F'
        gamma = 40.078;
    case '31P'
        gamma = 17.235;
    otherwise
        error('Nucleus %s not supported yet.', nucleus);
end

end