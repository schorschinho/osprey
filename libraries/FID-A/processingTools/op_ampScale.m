% op_ampScale.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_ampScale(in,A);
% 
% DESCRIPTION:
% Scale the amplitude of a spectrum by factor A.
% 
% INPUTS:
% in    = input data in matlab structure format
% A     = Amplitude scaling factor.
%
% OUTPUTS:
% out   = Output following amplitude scaling.  

function out=op_ampScale(in,A);

sz_A = length(A);
out=in;
if sz_A == 1
    out.specs=in.specs*A;
    out.fids=in.fids*A;
else
    
    if in.dims.averages==0   
        for ss = 1 : out.sz(end)
            out.specs(:,ss)=in.specs(:,ss)*A(ss);
            out.fids(:,ss)=in.fids(:,ss)*A(ss);
        end
    else
        for ss = 1 : out.sz(end)
            out.specs(:,:,ss)=in.specs(:,:,ss)*A(ss);
            out.fids(:,:,ss)=in.fids(:,:,ss)*A(ss);
        end
    end
end

