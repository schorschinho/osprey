function [out,switchOrder] = osp_onOffClassifyMEGA(in, target)
%% [outA,switchOrder] = osp_onOffClassifyMEGA(in,target)
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
%       [out] = osp_onOffClassifyMEGA(in target)
%
%   INPUTS:
%       in     = FID-A structure containing MEGA spectrum.
%       target  = String. Can be 'GABA' or 'GSH'.
%
%   OUTPUTS:
%       out    = FID-A structure containing the ordered MEGA sub-spectrum.
%       switchOrder = Vector indicating the order in which the input
%               spectra are rearranged in order to generate the output order.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-03-18)
%       goeltzs1@jhmi.edu
%           
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-03-18: First version of the code.
%       2020-01-15: Modified by Helge Zoellner (absolute value
%       implementation)


switch target
    case 'GABA'
        % Determine which of the differences has an upright NAA peak
        temp = op_freqrange(in, 1.7, 2.3);        
    case 'GSH'
        % Determine which of the differences has an upright water peak
        temp = op_freqrange(in, 4, 5);
    case 'Lac'
        % Determine which of the differences has an upright water peak
        temp = op_freqrange(in, 3.8, 4);
    case 'PE398'
        % Determine which of the differences has an upright water peak
        temp = op_freqrange(in, 3.8, 4.1);       
    case 'PE322'
        % Determine which of the differences has an upright water peak
        temp = op_freqrange(in, 3.0, 3.3);     
    otherwise
        disp('MEGA ON/OFF classifier does not recognize the input argument ''target''. We automatically assume no reordering. You can change that in osp_onOFFClassifyMEGA.m');
         outA = inA;
        outB = inB;
        switchOrder = 0;
end

spec = abs(real(temp.specs));        
max_diffAB = max(spec(:,1) - spec(:,2));
max_diffBA = max(spec(:,2) - spec(:,1));
out = in;

if max_diffAB > max_diffBA

    switchOrder = 0;
else
    temp = in;
    out.specs(:,1) = temp.specs(:,2);
    out.specs(:,2) = temp.specs(:,1);
    out.fids(:,1) = temp.fids(:,2);
    out.fids(:,2) = temp.fids(:,1);
    switchOrder = 1;
end

end