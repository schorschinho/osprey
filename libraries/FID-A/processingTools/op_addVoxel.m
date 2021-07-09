% op_addVoxel.m
% Helge Zoellner, The Johns Hopkins University 2020.
% 
% USAGE:
% out=op_takeVoxel(in1,in2,index);
% 
% DESCRIPTION:
% Adds  a certain voxel with indices corresponding to the 'index' input
% array. For multivoxel it is a single input, for 2D-MRSI a vector with 2 indices 
% corresponding to the x- and y-directions, and for 3D-MRSI a vector with 3 indices
% corresponding to the x-,y-, and z-directions.
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the subspectra you would like 
%          to extract.
%
% OUTPUTS:
% out    = Output dataset consisting of voxel index extracted from 
%          the input.

function out=op_addVoxel(in1,in2,index);

fids=in1.fids(:,:,:,:,:,:,:,:,:);

%change the dims variables
dims = in1.dims;
if length(index)==1 
    if ~isfield(dims, 'Xvoxels')
        dims.Xvoxels=length(size(fids));
    else if dims.Xvoxels == 0
            dims.Xvoxels=length(size(fids));
        end
    end
end
if length(index)==2 
    if ~isfield(dims, 'Xvoxels')
        dims.Xvoxels=length(size(fids));
        dims.Yvoxels=length(size(fids))+1;    
    else if dims.Xvoxels == 0
            dims.Xvoxels=length(size(fids));
            dims.Yvoxels=length(size(fids))+1; 
        end
    end
end
if length(index)==3 
    if ~isfield(dims, 'Xvoxels')
        dims.Xvoxels=length(size(fids));
        dims.Yvoxels=length(size(fids))+1;
        dims.Zvoxels=length(size(fids))+2;    
    else if dims.Xvoxels == 0
            dims.Xvoxels=length(size(fids));
            dims.Yvoxels=length(size(fids))+1;
            dims.Zvoxels=length(size(fids))+2; 
        end
    end        
end

%PRIAM data
if length(index)==1
    if dims.Xvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif dims.Xvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1.  This should never happen!  Aborting!');
    elseif dims.Xvoxels==2
        fids(:,index)=in2.fids(:);
    elseif dims.Xvoxels==3;
        fids(:,:,index)=in2.fids(:,:);
    elseif dims.Xvoxels==4;
        fids(:,:,:,index)=in2.fids(:,:,:);
    elseif dims.Xvoxels==5
        fids(:,:,:,:,index)=in2.fids(:,:,:,:);
    elseif dims.Xvoxels==6
        fids(:,:,:,:,:,index)=in2.fids(:,:,:,:,:);    
    end
end

% 2D MRSI data
if length(index)==2
    if dims.Xvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif dims.Xvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1 or dims.Yvoxels==1.  This should never happen!  Aborting!');
    elseif dims.Xvoxels==2
        fids(:,index(1),index(2))=in2.fids(:);
    elseif dims.Xvoxels==3
        fids(:,:,index(1),index(2))=in2.fids(:,:);
    elseif dims.Xvoxels==4
        fids(:,:,:,index(1),index(2))=in2.fids(:,:,:);
    elseif dims.Xvoxels==5
        fids(:,:,:,:,index(1),index(2))=in2.fids(:,:,:,:);
    elseif dims.Xvoxels==6
        fids(:,:,:,:,:,index(1),index(2))=in2.fids(:,:,:,:,:);    
    end
end

% 3D MRSI data
if length(index)==3
    if dims.Xvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif dims.Xvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1 or dims.Yvoxels==1 or dims.Zvoxels==1.  This should never happen!  Aborting!');
    elseif dims.Xvoxels==2
        fids(:,index(1),index(2),index(3))=in2.fids(:);
    elseif dims.Xvoxels==3
        fids(:,:,index(1),index(2),index(3))=in2.fids(:,:);
    elseif dims.Xvoxels==4
        fids(:,:,:,index(1),index(2),index(3))=in2.fids(:,:,:);
    elseif dims.Xvoxels==5
        fids(:,:,:,:,index(1),index(2),index(3))=in2.fids(:,:,:,:);
    elseif dims.Xvoxels==6
        fids(:,:,:,:,:,index(1),index(2),index(3))=in2.fids(:,:,:,:,:);  
    end
end

%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in1.dims.t),in1.dims.t);


%re-calculate the sz variable
sz=size(fids);


%FILLING IN DATA STRUCTURE
out=in1;
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

%Adding processed specific fields
fields = {'specReg','target','watersupp'};
for f = 1 : length(fields)
    if isfield(out,fields{f})
        if iscell(out.(fields{f}))
            %PRIAM data
            if length(index)==1
                out.(fields{f}){index} = in2.(fields{f});
            end

            % 2D MRSI data
            if length(index)==2
                out.(fields{f}){index(1),index(2)} = in2.(fields{f});
            end

            % 3D MRSI data
            if length(index)==3
                out.(fields{f}){index(1),index(2),index(3)} = in2.(fields{f});
            end
        else
            if length(index)==1
                out = rmfield(out, fields{f});
                out.(fields{f}){1} = in1.(fields{f});
                out.(fields{f}){index} = in2.(fields{f});
            end

            % 2D MRSI data
            if length(index)==2
                out = rmfield(out, fields{f});
                out.(fields{f}){1,1} = in1.(fields{f});
                out.(fields{f}){index(1),index(2)} = in2.(fields{f});
            end

            % 3D MRSI data
            if length(index)==3
                out = rmfield(out, fields{f});
                out.(fields{f}){1,1,1} = in1.(fields{f});
                out.(fields{f}){index(1),index(2),index(3)} = in2.(fields{f});
            end
        end            
    end
end

%FILLING IN THE FLAGS
out.flags=in1.flags;
out.flags.MultiVoxel=1;