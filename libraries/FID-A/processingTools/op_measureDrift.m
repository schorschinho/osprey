% op_measureDrift.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [out] = op_measureDrift(in);
%
% DESCRIPTION:
% Measures the frequency drift of the Cr signal over the course of the
% acquisition.
%
% INPUTS:
% in        = Input data structure.
%
% OUTPUTS:
% out       = Vector containing the frame-by-frame maxima of the Cr signal.


function [out] = op_measureDrift(in)

% Check whether data is coil-combined. If not, throw error.
if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% Check whether data is averaged. If it is, return zero.
if in.dims.averages == 0 || in.averages < 2
    out = 0;
    return
end

% Check whether input structure has sub-spectra.
% If it does, prepare to loop over all sub-spectra. The drift will be
% measured for each subspectrum separately.
if in.dims.subSpecs == 0
    numSubSpecs = 1;
else
    numSubSpecs = in.sz(in.dims.subSpecs);
end

% Pre-allocate memory.
nAvgs   = in.sz(in.dims.averages);
out     = zeros(nAvgs,numSubSpecs);

for mm=1:numSubSpecs
    % Phase the input spectrum by fitting a Cr-Cho doublet
    [in,~]       = op_phaseCrCho(in, 1);
    
    % For each average, fit a Cr-Cho doublet
    for kk = 1:nAvgs
        temp_kk = op_takeaverages(in, kk);
        parsFit = op_creChoFit(temp_kk, 1);
        out(kk,mm) = parsFit(3);
    end
end

end