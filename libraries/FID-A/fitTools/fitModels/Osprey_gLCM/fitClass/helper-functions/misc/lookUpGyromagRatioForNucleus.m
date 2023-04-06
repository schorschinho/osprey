function gamma = lookUpGyromagRatioForNucleus(nucleus)
%% gamma = lookUpGyromagRatioForNucleus(nucleus)
% Returns the gyromagnetic ratio [MHz/T] depending
% on the nucleus.
%
%   USAGE:
%       gamma = lookUpGyromagRatioForNucleus(nucleus)
%
%   INPUTS:
%       in             =  string with nucleus name
%
%   OUTPUTS:
%       gamma          = gyromagnetic ratio [MHz/T]
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

if iscell(nucleus)                  % Is cell
    nucleus = nucleus{1};           % Get first entry
end

switch strtrim(nucleus)             % Switch nucleus string
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