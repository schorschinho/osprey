% op_HadamardScans.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_addScans(in1,combMatrix);
% 
% DESCRIPTION:
% Add or subtract two scans together.
% 
% INPUTS:
% in        = Spectrum (in matlab structure format)
% combMatrix	 = How to combine the spectra
%
% OUTPUTS:
% out        = Result of adding inputs in1 and in2.  


function out=op_HadamardScans(in,combMatrix,name);

dim = in.dims.subSpecs;

if ~isempty(fieldnames(in))
    if dim == 2
        if length(combMatrix) == 2
            fids = combMatrix(1)*in.fids(:,1) + combMatrix(2)*in.fids(:,2);
        end
        if length(combMatrix) == 4
            fids = combMatrix(1)*in.fids(:,1) + combMatrix(2)*in.fids(:,2) + combMatrix(3)*in.fids(:,3)+ combMatrix(4)*in.fids(:,4);
        end
    end
end

%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

out=in;
if ~contains(in.names,name)
        out.fids=cat(dim,in.fids,fids);
        out.specs=cat(dim,in.specs,specs);
        out.sz=size(fids);            
        out.names{end+1} = name;
        out.subspecs = out.subspecs + 1;
    else
        ind = find(strcmp(out.names,name));
        out.fids(:,ind,:,:,:,:,:) = fids(:,1,:,:,:,:,:);
        out.specs(:,ind,:,:,:,:,:) = specs(:,1,:,:,:,:,:);
end
       



