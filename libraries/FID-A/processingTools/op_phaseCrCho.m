% op_phaseCrCho.m
% Georg Oeltzschner, Johns Hopkins University 2019.
% 
% USAGE:
% [out,phaseShift]=op_phaseCrCho(in, suppressPlot);
% 
% DESCRIPTION:
% This function fits two Lorentzians to the Cr and Cho peaks, and optimizes
% the phase of the input spectrum based on Cr/Cho.
% 
% INPUTS:
% in            = input data in matlab structure format. 
% suppressPlot  = (optional) Boolian to suppress plots.  Default = 0;
%
% OUTPUTS:
% out        = Output following automatic phasing.
% phShift    = The phase shift (in degrees) that was applied.

function [out,phShift]=op_phaseCrCho(in, suppressPlot);

if nargin < 2
    suppressPlot = 0;
end

if in.dims.coils>0
    error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
end
if in.dims.averages>0
    error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
end
if in.dims.extras>0
    error('ERROR:  Can not operate on data with extras dimension!  ABORTING!!');
end
if in.dims.subSpecs>0
    error('ERROR:  Can not operate on data with multiple subspectra!  ABORTING!!');
end

% Fit the Cr-Cho doublet to the input spectrum
parsFit = op_creChoFit(in, suppressPlot);
% Extract phase adjustment
phShift = parsFit(4);
% Apply negative phase of the fit to the data
out     = op_addphase(in, -parsFit(4)*180/pi, 0, in.centerFreq, suppressPlot);

end

