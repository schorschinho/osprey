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
% If data with multiple averages is supplied, the fit is performed to the
% average, and all averages are phased accordingly.
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
    in_avg = op_averaging(in);
else
    in_avg = in;
end
if in.dims.extras>0
    error('ERROR:  Can not operate on data with extras dimension!  ABORTING!!');
end
if in.dims.subSpecs>0
    temp = cell(in.subspecs,1);
    phShift = zeros(in.subspecs,1);
    temp{1} = op_takesubspec(in_avg,1);
    % Fit the Cr-Cho doublet to the input spectrum
    parsFit = op_creChoFit(temp{1}, suppressPlot);
    % Extract phase adjustment
    phShift = ones(in.subspecs,1) * parsFit(4);
        % Apply negative phase of the fit to the data
    for ss = 1 : in.subspecs   
        temp{ss} = op_takesubspec(in_avg,ss);
        temp{ss}     = op_addphase(temp{ss}, -phShift(ss)*180/pi, 0, temp{ss}.centerFreq, suppressPlot);   
    end
    out = temp{1};
    for ss = 2 : in.subspecs
        out = op_mergesubspec(out,temp{ss});
    end
    
else
    % Fit the Cr-Cho doublet to the input spectrum
    parsFit = op_creChoFit(in_avg, suppressPlot);
    % Extract phase adjustment
    phShift = parsFit(4);
    % Apply negative phase of the fit to the data
    out     = op_addphase(in, -parsFit(4)*180/pi, 0, in.centerFreq, suppressPlot);
   
end

% Add NIfTI-MRS provenance
% Generate fields for provenance
fields.Method   = 'Phasing';
fields.Details  = ['Phase spectrum using a two Lorentzian Cr/Cho fit'];
in = op_add_analysis_provenance(in,fields);
end

