% op_sortMRSIsdat.m
% Helge Zoellner, The Johns Hopkins University 2020.
% 
% USAGE:
% out=osp_MRSIRecon(in);
% 
% DESCRIPTION:
%
% 
% in         = input data in matlab structure format.
%
% OUTPUTS:
% out        = Output dataset following zeropadding.

function out=op_sortMRSIsdat(in);

% 2D MRSI data
%Now that we've indexed the dimensions of the data array, as if they were
%single voxel data. This means taht the different voxels are stored in the
% averages dimensions. We need to re-arange it:  
%   1) time domain data.  
%   2) coils.
%   3) averages.
%   4) subSpecs.
%   5) extras.
%   6) Xvoxels
%   7) Yvoxels

% Create a number code for the switch 
dims = in.dims;
if ~isfield(dims, 'Zvoxels')
    dimsCode = [in.dims.t,in.dims.coils,in.dims.averages,in.dims.subSpecs,in.dims.extras];
else
    dimsCode = [in.dims.t,in.dims.coils,in.dims.averages,in.dims.subSpecs,in.dims.extras,in.dims.Zvoxels];
end
dimsCode(dimsCode>0) = 1;
dimsCode = num2str(dimsCode);
dimsCode(dimsCode == ' ') = [];

%Calculate averages from unsorted data
if in.nZvoxels <= 1
    averages = in.sz(dims.averages)/(in.nXvoxels * in.nYvoxels);
else if isfield(dims, 'Zvoxels')
    averages = in.sz(dims.averages)/(in.nXvoxels * in.nYvoxels);
    else
    averages = in.sz(dims.averages)/(in.nXvoxels * in.nYvoxels* in.nZvoxels);
    end
end

if in.nZvoxels <= 1
    % 2D MRSI data
    switch dimsCode
        case '10100'
            fids = reshape(in.fids,[in.sz(dims.t) averages in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=3;
            dims.Yvoxels=4;  
        case '10101'
            fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=4;
            dims.Yvoxels=5; 
        case '10110'
            fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=4;
            dims.Yvoxels=5; 
        case '10111'
            fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=5;
            dims.Yvoxels=6; 
        case '11101'
            fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=6;
            dims.Yvoxels=7; 
        case '11110'
            fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=6;
            dims.Yvoxels=7; 
        case '11111'
            fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels]);
            dims.Xvoxels=8;
            dims.Yvoxels=9; 
    end
end
% 3D MRSI data
if in.nZvoxels > 1
    if ~isfield(dims, 'Zvoxels')
        switch dimsCode
            case '10100'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.nXvoxels in.nYvoxels in.nZvoxels]);
                dims.Xvoxels=4;
                dims.Yvoxels=5;
                dims.Zvoxels=6; 
            case '10101'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=4;
                dims.Yvoxels=5;
                dims.Zvoxels=6; 
            case '10110'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels in.nZvoxels]);            
                dims.Xvoxels=4;
                dims.Yvoxels=5;
                dims.Zvoxels=6; 
            case '10111'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '11101'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '11110'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '11111'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=6;
                dims.Yvoxels=7;
                dims.Zvoxels=8; 
        end
    else
        switch dimsCode
            case '101001'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.nXvoxels in.nYvoxels in.nZvoxels]);
                dims.Xvoxels=3;
                dims.Yvoxels=4;
                dims.Zvoxels=5; 
            case '101011'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=4;
                dims.Yvoxels=5;
                dims.Zvoxels=6; 
            case '101101'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels in.nZvoxels]);            
                dims.Xvoxels=4;
                dims.Yvoxels=5;
                dims.Zvoxels=6; 
            case '101111'
                fids = reshape(in.fids,[in.sz(dims.t) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '111011'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '111101'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=5;
                dims.Yvoxels=6;
                dims.Zvoxels=7; 
            case '111111'
                fids = reshape(in.fids,[in.sz(dims.t) in.sz(dims.coils) averages in.sz(dims.subSpecs) in.sz(dims.extras) in.nXvoxels in.nYvoxels in.nZvoxels]);           
                dims.Xvoxels=6;
                dims.Yvoxels=7;
                dims.Zvoxels=8; 
        end
    end
end


%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);


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
out.flags.MultiVoxel=1;
end