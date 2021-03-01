function [varargout] = osp_editSubSpecAlignLNorm(varargin)
%% [varargout] = osp_editSubSpecAlign(varargin)
%   Aligns sub-spectra of edited MRS data to minimize
%   subtraction artefacts. It uses a L1 Norm minimization including the
%   whole frequency range (Cleve et al., JMR, 2017 (https://doi.org/10.1016/j.jmr.2017.04.004) .
%
%   USAGE:
%       [outA, outB, outC, outD] = osp_editSubSpecAlign(inA, inB, inC, inD);
%
%   INPUTS:
%       inA        = Input data structure with sub-spectrum A.
%       inB        = Input data structure with sub-spectrum B.
%       inC        = Input data structure with sub-spectrum C. (optional)
%       inD        = Input data structure with sub-spectrum D. (optional)
%
%   OUTPUTS:
%       outA       = Output following alignment of averages.
%       outB       = Output following alignment of averages.
%       outC       = Output following alignment of averages. (optional)
%       outD       = Output following alignment of averages. (optional)
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2021-03-01)
%       hzoelln2@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-01-01: First version of the code.

% Determine whether there are 2 (MEGA) or 4 (HERMES/HERCULES) inputs
if nargin == 4
    seqType = 'HERMES';
    inA     = varargin{1};
    inB     = varargin{2};
    inC     = varargin{3};
    inD     = varargin{4};
elseif nargin == 2
    seqType = 'MEGA';
    inA     = varargin{1};
    inB     = varargin{2};
else
    error('Error in osp_editSubSpecAlign! Needs to have either 2 (for MEGA) or 4 (for HERMES/HERCULES) sub-spectra provided');
end

% Check whether data is coil-combined. If not, throw error.
if ~inA.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% Check whether data is averaged. If not, throw error.
if ~inA.flags.averaged
    error('ERROR:  I think it only makes sense to do this after averaging using op_averaging.  ABORTING!!');
end

%%% 1. SET UP REQUIRED VARIABLES %%%
% Define the frequency ranges over which water, NAA, and Cho subtraction artefacts
% are to be minimized. Also, get a good starting estimate for the frequency
% alignment shift by determining the difference between the two maxima in
% the respective peak range (water, NAA, or Cho).
freq = inA.ppm;
freqLim(1,:) = freq <= 3.8 & freq >= 1.95;
x0(1,:) = [0 0];


% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Display','off','Algorithm','levenberg-marquardt');

% Initialize common variables
t           = inA.t;


%%% 2. PERFORM ALIGNMENT BASED ON SEQUENCE TYPE
if strcmp(seqType, 'HERMES')


    a = max(max([abs(real(inA.specs)) abs(real(inB.specs))]));
    fun = @(x) objFunc(op_ampScale(inA, 1/a), op_ampScale(inB, 1/a), freqLim(1,:), t, x);
    param(1,:) = lsqnonlin(fun, x0(1,:), [], [], lsqnonlinopts);
    % Apply the calculated frequency/phase adjustment to the inB spectrum
    fidsB = inB.fids.*exp(1i*pi*(t'*param(1,1)*2+param(1,2)/180));
    specsB = fftshift(fft(fidsB, [], inB.dims.t), inB.dims.t);
    % Create output
    outA = inA;
    outB = inB;
    outB.fids = fidsB;
    outB.specs = specsB;

    a = max(max([abs(real(inA.specs)) abs(real(inC.specs))]));
    fun = @(x) objFunc(op_ampScale(inA, 1/a), op_ampScale(inC, 1/a), freqLim(1,:), t, x);
    param(2,:) = lsqnonlin(fun, x0(1,:), [], [], lsqnonlinopts);
    % Apply the calculated frequency/phase adjustment to the inC spectrum
    fidsC = inC.fids.*exp(1i*pi*(t'*param(2,1)*2+param(2,2)/180));
    specsC = fftshift(fft(fidsC, [], inC.dims.t), inC.dims.t);
    % Create output
    outC = inC;
    outC.fids = fidsC;
    outC.specs = specsC;

    a = max(max([abs(real(outC.specs)) abs(real(inD.specs))]));
    fun = @(x) objFunc(op_ampScale(outC, 1/a), op_ampScale(inD, 1/a), freqLim(1,:), t, x);
    param(3,:) = lsqnonlin(fun, x0(1,:), [], [], lsqnonlinopts);
    % Apply the calculated frequency/phase adjustment to the inD spectrum
    fidsD = inD.fids.*exp(1i*pi*(t'*param(3,1)*2+param(3,2)/180));
    specsD = fftshift(fft(fidsD, [], inD.dims.t), inD.dims.t);

    % Create output
    outD = inD;
    outD.fids = fidsD;
    outD.specs = specsD;
        
         
    varargout{1} = outA;
    varargout{2} = outB;
    varargout{3} = outC;
    varargout{4} = outD;
    
elseif strcmp(seqType, 'MEGA')
    % For MEGA-edited data, the 'reporter signal' that is used to align the
    % two sub-spectra will depend on the edited metabolite. Since the peak
    % needs to be identical in both acquisitions, we choose the residual
    % water peak for GABA-edited data, and the NAA peak for GSH-edited
    % data.

    a = max(max([abs(real(inA.specs)) abs(real(inB.specs))]));
    fun = @(x) objFunc(op_ampScale(inA, 1/a), op_ampScale(inB, 1/a), freqLim(1,:), t, x);
    param(1,:) = lsqnonlin(fun, x0(1,:), [], [], lsqnonlinopts);

    % Apply the calculated frequency/phase adjustment to the inB spectrum
    fidsB = inB.fids.*exp(1i*pi*(t'*param(1,1)*2+param(1,2)/180));
    specsB = fftshift(fft(fidsB, [], inB.dims.t), inB.dims.t);
    % Create output
    outA = inA;
    outB = inB;
    outB.fids = fidsB;
    outB.specs = specsB;
    varargout{1} = outA;
    varargout{2} = outB;
    
end

end

function out = objFunc(in1, in2, freqLim, t, x)

% This is the objective function that minimizes the absolute sum of the difference over the 
% frequency range (freqLim) between the target spectrum (in1) and the 
% spectrum that is to be frequency-and-phase aligned (in2).

f   = x(1); % frequency correction
phi = x(2); % phase correction

y1 = in1.fids;
y2 = in2.fids .* exp(1i*pi*(t'*f*2+phi/180)); % apply to time-domain data

% fft
a = real(fftshift(fft(y1,[],1),1));
b = real(fftshift(fft(y2,[],1),1));

% return difference vector over defined frequency range
DIFF = a - b;
out = sum(abs(DIFF(freqLim)));

end