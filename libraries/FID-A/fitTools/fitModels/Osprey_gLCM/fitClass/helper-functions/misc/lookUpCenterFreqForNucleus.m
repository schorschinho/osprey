function centerFreq = lookUpCenterFreqForNucleus(nucleus)
% centerFreq = lookUpCenterFreqForNucleus(nucleus)
% Returns the ppm value at the center of the spectrum depending
% on the nucleus.
%
%   USAGE:
%       centerFreq = lookUpCenterFreqForNucleus(nucleus)
%
%   INPUTS:
%       in             =  string with nucleus name
%
%   OUTPUTS:
%       centerFreq     = center frequency [ppm]
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu

%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017) 
%%  Find gyromagnetic ratio

if iscell(nucleus)              % Is cell
    nucleus = nucleus{1};       % Get first entry
end

switch strtrim(nucleus)         % Switch nucleus string
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