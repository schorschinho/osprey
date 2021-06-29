% op_takeVoxel.m
% Helge Zoellner, The Johns Hopkins University 2020.
% 
% USAGE:
% out=op_takeVoxel(in,index);
% 
% DESCRIPTION:
% Extract the a certain voxel with indices corresponding to the 'index' input
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

function out=op_takeVoxel(in,index);



% if (in.flags.MultiVoxel==0)
%     error('ERROR:  This is not a multi voxel dataset!  Aborting!');
% end
if ~isfield(in.dims,'Zvoxels')
    index = index(1:2);
end
if ~isfield(in.dims,'Yvoxels')
    index = index(1);
end
%PRIAM data
if length(index)==1
    if in.dims.Xvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif in.dims.Xvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1.  This should never happen!  Aborting!');
    elseif in.dims.Xvoxels==2
        fids=in.fids(:,index);
    elseif in.dims.Xvoxels==3;
        fids=in.fids(:,:,index);
    elseif in.dims.Xvoxels==4;
        fids=in.fids(:,:,:,index);
    elseif in.dims.Xvoxels==5
        fids=in.fids(:,:,:,:,index);
        elseif in.dims.Xvoxels==6
        fids=in.fids(:,:,:,:,:,index);    
    end
end

% 2D MRSI data
if length(index)==2
    if in.dims.Xvoxels==0 || in.dims.Yvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif in.dims.Xvoxels==1 || in.dims.Yvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1 or dims.Yvoxels==1.  This should never happen!  Aborting!');
    elseif in.dims.Xvoxels==2
        fids=in.fids(:,index(1),index(2));
    elseif in.dims.Xvoxels==3;
        fids=in.fids(:,:,index(1),index(2));
    elseif in.dims.Xvoxels==4;
        fids=in.fids(:,:,:,index(1),index(2));
    elseif in.dims.Xvoxels==5
        fids=in.fids(:,:,:,:,index(1),index(2));
        elseif in.dims.Xvoxels==6
        fids=in.fids(:,:,:,:,:,index(1),index(2));    
    end
end

% 3D MRSI data
if length(index)==3
    if in.dims.Xvoxels==0 || in.dims.Yvoxels==0 || in.dims.Zvoxels==0
        %Can't take subspec because there are none:
        error('ERROR:  There are multiple voxels in this dataset!  Aborting!');
    elseif in.dims.Xvoxels==1 || in.dims.Yvoxels==1 || in.dims.Zvoxels==1
        %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
        error('ERROR:  dims.Xvoxels==1 or dims.Yvoxels==1 or dims.Zvoxels==1.  This should never happen!  Aborting!');
    elseif in.dims.Xvoxels==2
        fids=in.fids(:,index(1),index(2),index(3));
    elseif in.dims.Xvoxels==3;
        fids=in.fids(:,:,index(1),index(2),index(3));
    elseif in.dims.Xvoxels==4;
        fids=in.fids(:,:,:,index(1),index(2),index(3));
    elseif in.dims.Xvoxels==5
        fids=in.fids(:,:,:,:,index(1),index(2),index(3));
        elseif in.dims.Xvoxels==6
        fids=in.fids(:,:,:,:,:,index(1),index(2),index(3));  
    end
end

%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims = in.dims;
if length(index)==1
    dims.Xvoxels=0;
end
if length(index)==2
    dims.Xvoxels=0;
    dims.Yvoxels=0;    
end
if length(index)==3
    dims.Xvoxels=0;
    dims.Yvoxels=0;
    dims.Zvoxels=0;    
end

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

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.MultiVoxel=0;
end