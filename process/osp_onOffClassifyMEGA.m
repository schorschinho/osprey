function [outA, outB, switchOrder] = osp_onOffClassifyMEGA(inA, inB, target)
%% [outA, outB, switchOrder] = osp_onOffClassifyMEGA(inA, inB, target)
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
        tempA = op_freqrange(inA, 1.7, 2.3);
        tempB = op_freqrange(inB, 1.7, 2.3);

        specA = abs(real(tempA.specs));
        specB = abs(real(tempB.specs));
        
        max_diffAB = max(specA - specB);
        max_diffBA = max(specB - specA);
        
        
        if max_diffAB > max_diffBA
            outA = inA;
            outB = inB;
            switchOrder = 0;
        else
            outA = inB;
            outB = inA;
            switchOrder = 1;
        end
    case 'GSH'
        % Determine which of the differences has an upright water peak
        tempA = op_freqrange(inA, 4, 5);
        tempB = op_freqrange(inB, 4, 5);

        specA = abs(tempA.specs);
        specB = abs(tempB.specs);
        
        max_diffAB = max(specA - specB);
        max_diffBA = max(specB - specA);
        
        
        if max_diffAB > max_diffBA
            outA = inA;
            outB = inB;
            switchOrder = 0;
        else
            outA = inB;
            outB = inA;
            switchOrder = 1;
        end
    case 'Lac'
        % Determine which of the differences has an upright water peak
        tempA = op_freqrange(inA, 3.8, 4);
        tempB = op_freqrange(inB, 3.8, 4);

        specA = abs(tempA.specs);
        specB = abs(tempB.specs);
        
        max_diffAB = max(specA - specB);
        max_diffBA = max(specB - specA);
        
        
        if max_diffAB < max_diffBA
            outA = inB;
            outB = inA;
            switchOrder = 1;
        else
            outA = inA;
            outB = inB;
            switchOrder = 0;
        end 
    case 'PE398'
        % Determine which of the differences has an upright water peak
        tempA = op_freqrange(inA, 3.8, 4.1);
        tempB = op_freqrange(inB, 3.8, 4.1);

        specA = abs(tempA.specs);
        specB = abs(tempB.specs);
        
        max_diffAB = max(specA - specB);
        max_diffBA = max(specB - specA);
        
        
        if max_diffAB < max_diffBA
            outA = inB;
            outB = inA;
            switchOrder = 1;
        else
            outA = inA;
            outB = inB;
            switchOrder = 0;
        end   
        
    case 'PE322'
        % Determine which of the differences has an upright water peak
        tempA = op_freqrange(inA, 3.0, 3.3);
        tempB = op_freqrange(inB, 3.0, 3.3);

        specA = abs(tempA.specs);
        specB = abs(tempB.specs);
        
        max_diffAB = max(specA - specB);
        max_diffBA = max(specB - specA);
        
        
        if max_diffAB < max_diffBA
            outA = inB;
            outB = inA;
            switchOrder = 1;
        else
            outA = inA;
            outB = inB;
            switchOrder = 0;
        end        
    otherwise
        error('MEGA ON/OFF classifier does not recognize the input argument ''target''. Set to ''GABA'' or ''GSH''.');
end


end