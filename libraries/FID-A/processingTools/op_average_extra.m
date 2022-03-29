% op_averaging.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_average_extra(in,which);
% 
% DESCRIPTION:
% Combine the averages in a scan by adding the averages together and then 
% dividing by the number of averages.
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Output following averaging.  

function out=op_average_extra(in,to_average);
if nargin < 2
    to_average = 1;
end
if in.dims.extras==0 || in.extras<2
    %DO NOTHING
    warning('Just one extra spectrum found. Returning input without modification!');
    out=in;
else
    if to_average == 1
        %add the spectrum along the averages dimension;
        fids=sum(in.fids,in.dims.extras);
        fids=squeeze(fids);
        fids=fids/in.sz(in.dims.extras); %divide by number of averages;

        %re-calculate Specs using fft
        specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

        %change the dims variables.
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
        if in.dims.subSpecs>in.dims.extras
            dims.subSpecs=in.dims.subSpecs-1;
        else
            dims.subSpecs=in.dims.subSpecs;
        end
        dims.extras = 0;        
        
    else
        
    end

    
    
    %re-calculate the sz variable
    sz=size(fids);
    
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.dims=dims;
    out.averages=1;
    
    %Change all the remaining entries
    if to_average == 1
        out.spectralwidth = mean(out.spectralwidth);
        out.dwelltime = mean(out.dwelltime);
        out.txfrq = mean(out.txfrq);
        out.seq = out.seq(1);
        out.te = mean(out.te);
        out.tr = mean(out.tr);
        out.pointsToLeftshift = mean(out.pointsToLeftshift);
        out.centerFreq = mean(out.centerFreq);
        out.names = out.names(1,:);
        out.extra_names = {'Average across all extra entries'};
        out.extras = 0;
    else
        
    end
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.averaged=1;

end



