% op_mergesubspec.m
% Helge Zoellner, Johns Hopkins University 2021.
% 
% USAGE:
% out=op_mergeextra(varargin);
% 
% DESCRIPTION:
% Merges input into the extra dimension
% 
% INPUTS:
% varargin     = input data in matlab structure format.
%
% OUTPUTS:
% out    = Output dataset with merged subspectra

function out=op_mergeextra(varargin);
in1 = varargin{1};

if isempty(varargin{2}) %Check if the first argument to be a null string
    error('You have to at least put in two structs')
else
    if isfield(in1, 'extra_names')
        names = in1.extra_names;
    else
        names = {};
    end
    if isfield(in1, 'exp_var')
        exp_var = in1.exp_var;
    else
        exp_var = 1;
    end
    if ~ischar(varargin{end})
        error('You have to at indicate an extra name as last argument')     
    else
        names{end+1} = varargin{end};
    end
    
    if ~isfield(in1, 'extras')
        in1.extras = 1;
    end
    in1.extra_names = names;
    in1.exp_var = exp_var;
    
    if ischar(in1.seq)
        temp_seq = in1.seq;
        in1 = rmfield(in1,'seq');
        in1.seq{1} = temp_seq;
    end
    for ex = 2 : length(varargin)-1
        in2 = varargin{ex};
        if in1.dims.subSpecs~=0 && in2.dims.subSpecs~=0 && in1.dims.subSpecs ~= in2.dims.subSpecs
            %RULE:  Averages dimension must be the same for both inputs UNLESS one of
            %the inputs has an averages dimension of zero:
            error('subspectra dimensions must be the same for both inputs');
        end

        if in1.dims.subSpecs == 2
            if in1.dims.extras==0 && in2.dims.extras==0
                dim=3;
            elseif in1.dims.extras==0
                dim=in2.dims.extras;
            elseif in2.dims.averages==0
                dim=in1.dims.extras;
            else
                dim=in1.dims.extras;
            end
        else
            if in1.dims.extras==0 && in2.dims.extras==0
                dim=2;
            elseif in1.dims.extras==0
                dim=in2.dims.extras;
            elseif in2.dims.averages==0
                dim=in1.dims.extras;
            else
                dim=in1.dims.extras;
            end
        end
        in1.dims.extras = dim;
        in1.fids=cat(dim,in1.fids,in2.fids);
        in1.specs=cat(dim,in1.specs,in2.specs);
        if isfield(in2,'extras') && (in2.extras > 0)
            in1.extras=in1.extras+in2.extras;
        else
            in1.extras=in1.extras+1;
        end
        in1.seq{end+1} = in2.seq;       
        in1.spectralwidth(end+1) = in2.spectralwidth;
        in1.dwelltime(end+1) = in2.dwelltime;
        in1.txfrq(end+1) = in2.txfrq;
        in1.te(end+1) = in2.te;
        in1.tr(end+1) = in2.tr;
        in1.centerFreq(end+1) = in2.centerFreq;
        if isfield(in1, 'pointsToLeftshift')
            in1.pointsToLeftshift(end+1) = in2.pointsToLeftshift;
        end
        if isfield(in1,'specReg')
            in1.specReg{end+1} = in2.specReg;
        end
        if isfield(in2,'exp_var')
            in1.exp_var(end+1) = in2.exp_var;
        else
            in1.exp_var(end+1) = 1;
        end
        if size(in1.names,2) > size(in2.names,2)
            in2.names = cat(2,in2.names,cell(1,size(in1.names,2)-size(in2.names,2)));
        else if size(in2.names,2) > size(in1.names,2)
            in1.names{size(in1.names,1),size(in2.names,2)} = [];
            end
        end
        if isfield(in1,'watersupp')
            if size(in1.watersupp,2) > size(in2.watersupp,2)
                in2.watersupp = cat(2,in2.watersupp,cell(1,size(in1.watersupp,2)-size(in2.watersupp,2)));
            else if size(in2.watersupp,2) > size(in1.watersupp,2)
                in1.watersupp{size(in1.watersupp,1),size(in2.watersupp,2)} = [];
                end
            end        
            in1.watersupp = cat(1,in1.watersupp,in2.watersupp);
        end
        in1.names = cat(1,in1.names,in2.names);
    end
    sz=size(in1.fids);
    %FILLING IN DATA STRUCTURE
    out=in1;
    out.sz=sz;

    %FILLING IN THE FLAGS
    out.flags=in1.flags;
    out.flags.writtentostruct=1;
end
