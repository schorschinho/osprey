% op_addextra.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_takesubspec(in,index);
% 
% DESCRIPTION:
% Add the subspectra with indices corresponding to the 'index' input
% array. 
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the subspectra you would like 
%          to extract.
%
% OUTPUTS:
% out    = Output dataset consisting of subspectra indices extracted from 
%          the input.

function out=op_addextra(in,in2,index);

if in.dims.extras>0
in.fids = squeeze(in.fids);        
    if in.dims.extras==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.extras==1.  This should never happen!  Aborting!');
    elseif in.dims.extras==2
        in.fids(:,index,:,:,:) = in2.fids(:,1,:,:,:);
    elseif in.dims.extras==3;
        in.fids(:,:,index,:,:,:) =  in2.fids(:,:,1,:,:,:);
    elseif in.dims.extras==4;
        in.fids(:,:,:,index,:,:,:) =  in2.fids(:,:,:,1,:,:,:);
    elseif in.dims.extras==5
        in.fids(:,:,:,:,index,:,:,:) =  in2.fids(:,:,:,:,1,:,:,:);
    end

    %re-calculate Specs using fft
    in.specs=fftshift(fft(in.fids,[],in.dims.t),in.dims.t);

    %re-calculate the sz variable
    sz=size(in.fids);

    %add manual phasing entries
    if isfield(in2, 'manual')
        in.manual = in2.manual;
    end

    %FILLING IN DATA STRUCTURE
    out=in;
    out.sz=sz;
    
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.isISIS=0;
else
    out = in;
end
