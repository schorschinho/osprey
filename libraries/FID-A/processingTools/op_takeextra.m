% op_takeextra.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_takesubspec(in,index);
% 
% DESCRIPTION:
% Extract the subspectra with indices corresponding to the 'index' input
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

function out=op_takeextra(in,index);

if in.dims.extras>0
in.fids = squeeze(in.fids);        
    if in.dims.extras==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.extras==1.  This should never happen!  Aborting!');
    elseif in.dims.extras==2
        fids=in.fids(:,index,:,:,:);
    elseif in.dims.extras==3;
        fids=in.fids(:,:,index,:,:,:);
    elseif in.dims.extras==4;
        fids=in.fids(:,:,:,index,:,:,:);
    elseif in.dims.extras==5
        fids=in.fids(:,:,:,:,index,:,:,:);
    end

    %re-calculate Specs using fft
    specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

    %change the dims variables
    if in.dims.t>in.dims.extras
        dims.t=in.dims.t-1;
    else
        dims.t=in.dims.t;
    end
    if in.dims.coils>in.dims.extras
        dims.coils=in.dims.coils-1;
    else
        dims.coils=in.dims.coils;
    end
    if in.dims.averages>in.dims.extras
        dims.averages=in.dims.averages-1;
    else
        dims.averages=in.dims.averages;
    end
    if in.dims.extras>in.dims.extras
        dims.subSpecs=in.dims.subSpecs-1;
    else
        dims.subSpecs=in.dims.subSpecs;
    end
    dims.extras=0;

    %re-calculate the sz variable
    sz=size(fids);


    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.dims=dims;
    if ~out.flags.averaged
        try
            out.averages = out.sz(out.dims.averages);
        catch
            out.averages = out.sz(2);
        end
    end
    if (out.dims.averages == 0) && (out.dims.extras == 0)
        out.dims.averages = 2;
    end
    out.extras=1;
    if isfield(out, 'extra_names')
        out.extra_names = out.extra_names(index);
    end
    if isfield(out, 'names')
        try
            out.names = out.names(index);
        catch
            out.names = out.names(1);
        end
    end
    if isfield(out, 'exp_var')
        out.exp_var = out.exp_var(index);
    end
    out.spectralwidth = out.spectralwidth(index);
    out.dwelltime = out.dwelltime(index);
    out.txfrq = out.txfrq(index);
    out.seq = out.seq{index};
    out.te = out.te(index);
    out.tr = out.tr(index);
    if isfield(out, 'pointsToLeftshift')
        out.pointsToLeftshift = out.pointsToLeftshift(index);
    end
    out.centerFreq = out.centerFreq(index);
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.isISIS=0;
else
    out = in;
end
