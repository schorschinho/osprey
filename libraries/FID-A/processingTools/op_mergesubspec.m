% op_mergesubspec.m
% Helge Zoellner, Johns Hopkins University 2021.
% 
% USAGE:
% out=op_mergesubspec(varargin);
% 
% DESCRIPTION:
% Merges input structs into the subspectra dimension
% 
% INPUTS:
% varargin     = input data in matlab structure format.
%
% OUTPUTS:
% out    = Output dataset with merged subspectra

function out=op_mergesubspec(varargin);
in1 = varargin{1};

if isempty(varargin{2}) %Check if the first argument to be a null string
    error('You have to at least put in two structs')
else
    
    for ss = 2 : length(varargin)
        in2 = varargin{ss};
        if in1.dims.subSpecs~=0 && in2.dims.subSpecs~=0 && in1.dims.subSpecs ~= in2.dims.subSpecs
            %RULE:  Averages dimension must be the same for both inputs UNLESS one of
            %the inputs has an averages dimension of zero:
            error('subspectra dimensions must be the same for both inputs');
        end

        if in1.dims.subSpecs==0 && in2.dims.subSpecs==0
            dim=2;
        elseif in1.dims.subSpecs==0
            dim=in2.dims.subSpecs;
        elseif in2.dims.averages==0;
            dim=in1.dims.subSpecs;
        else
            dim=in1.dims.subSpecs;
        end
        if ~contains(in1.names,in2.names(1))
            in1.fids=cat(dim,in1.fids,in2.fids);
            in1.specs=cat(dim,in1.specs,in2.specs);
            in1.subspecs=in1.subspecs+in2.subspecs;
            in1.names{end+1} = in2.names{1};
        else
            ind = find(strcmp(in1.names,in2.names(1)));
            in1.fids(:,ind,:,:,:,:,:) = in2.fids(:,1,:,:,:,:,:);
            in1.specs(:,ind,:,:,:,:,:) = in2.specs(:,1,:,:,:,:,:);
        end
        
    end
    sz=size(in1.fids);
    
    %FILLING IN DATA STRUCTURE
    out=in1;
    out.sz=sz;
%     out.rawAverages=in1.rawAverages+in2.rawAverages;
    if in1.dims.subSpecs==0
        out.dims.subSpecs=dim;
    end
    
    
    %FILLING IN THE FLAGS
    out.flags=in1.flags;
    out.flags.writtentostruct=1;
end
