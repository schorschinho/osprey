%io_writesdat.m
%Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% RF=io_writesdat(in,outfile);
% 
% DESCRIPTION:
% Takes MRS data in matlab structure format and writes it to a SDAT file
% that can be read by LCModel.
% 
% INPUTS:
% in         = input data in matlab structure format.
% sdatFile   = Desired filename of output SDAT file.
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in SDAT format. 

function RF=io_writesdat(in,sdatFile);
%function RF=writesdat(in,sdatFile);

% if in.flags.isISIS
%     error('ERROR:  Must make subspecs first');
% end

if ~in.flags.averaged
    disp('WARNING:  Signals must be averaged first');
end

if ~in.flags.addedrcvrs
    error('ERROR:  receiver channels must be combined first');
end

% Create dummy copy to return
RF = in;

% Create dummy array to hold FIDs
fids = in.fids;
reshape_data = reshape(fids, [1 in.sz(1) in.sz(2)]);
data=[real(reshape_data);imag(reshape_data)];

% Open a file with the designated name, and write the FID
fid     = fopen(sdatFile,'w','ieee-le');
status  = fwriteVAXG(fid,data,'float32');
fclose(fid);