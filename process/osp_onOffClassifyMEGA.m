function [outA, outB] = osp_onOffClassifyMEGA(inA, inB, target)
%% [outA, outB] = osp_onOffClassifyMEGA(inA, inB, target)
%   This function decides which of the two provided MEGA sub-spectra in the
%   are the edit-ON or the edit-OFF.
%
%   To this end, the difference spectrum is analyzed. For GABA-edited data,
%   the NAA peak is investigated; for GSH-edited data, the residual water
%   peak.
%
%   The function then assigns the edit-OFF spectrum to field A, and the
%   edit-ON spectrum to field B.
%
%   USAGE:
%       [outA, outB] = osp_onOffClassifyMEGA(inA, inB, target)
%
%   INPUTS:
%       inA     = FID-A structure containing one MEGA sub-spectrum.
%       inB     = FID-A structure containing the other MEGA sub-spectrum.
%       target  = String. Can be 'GABA' or 'GSH'.
%
%   OUTPUTS:
%       outA    = FID-A structure containing the edit-OFF MEGA sub-spectrum.
%       outB    = FID-A structure containing the edit-ON MEGA sub-spectrum.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-03-18)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-03-18: First version of the code.

% Calculate difference spectrum
DIFF = op_addScans(inA, inB, 1);

switch target
    case 'GABA'
        % Determine whether DIFF spectrum has a positive NAA peak.
        parsFit = op_naaFit(DIFF,0,0);

        if (parsFit(1) * parsFit(2)) > 0 % this is the area under the curve of the NAA fit
            outA = inA;
            outB = inB;
        else
            outA = inB;
            outB = inA;
        end
    case 'GSH'
        error('Automatic ON/OFF classification for GSH-edited data coming soon!');
    otherwise
        error('MEGA ON/OFF classifier does not recognize the input argument ''target''. Set to ''GABA'' or ''GSH''.');
end

%% Clean up

% Close any remaining open figures
close all;

end
