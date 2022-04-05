% io_writeniimrs.m
% Georg Oeltzschner, Johns Hopkins University 2021.
%
% USAGE:
% nii = io_writeniimrs(in, outfile);
% 
% DESCRIPTION:
% Takes MRS data in matlab structure format and writes it to a NIfTI-MRS
% file.
% See the specification under
% https://docs.google.com/document/d/1tC4ugzGUPLoqHRGrWvOcGCuCh_Dogx_uu0cxKub0EsM/edit
%
% This function is currently work-in-progress and is being tested on more
% and more datasets. Please contact the FID-A developers for help with a
% particular dataset that you might encounter problems with when using
% this function.
%
% Currently, this function is limited to single-voxel MRS data. It is
% planned to develop it for full compatibility with 2D and 3D multi-voxel
% and spectroscopic imaging data.
%
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output NIfTI-MRS (*.nii) file.
% dim_names  = Desired output dimension names (dim_5_name, dim_6_name,
%              dim_7_name)
%
% DEPENDENCIES:
% This function requires the dcm2nii toolbox (Xiangrui Li) to be on the
% MATLAB path
% https://github.com/xiangruili/dicm2nii
%
% OUTPUTS:
% nii         = Same as input. Not used. 

function nii = io_writeniimrs(in, outfile, dim_names)
%function nii = io_writeniimrs(in, outfile)

% Bring FID array into NIfTI-MRS shape
% This will need work as we move to 2D and 3D arrays
dims = in.dims;
% Find non-singleton dimensions
sqzDims = {};
dimsFieldNames = fieldnames(dims);
for rr = 1:length(dimsFieldNames)
    if dims.(dimsFieldNames{rr}) ~= 0
        sqzDims{end+1} = dimsFieldNames{rr};
    end
end
    
if ~isfield(dims, 'x')
    dims.x = 0;
    dim_1  = 1;
end
if ~isfield(dims, 'y')
    dims.y = 0;
    dim_2  = 1;
end
if ~isfield(dims, 'z')
    dims.z = 0;
    dim_3  = 1;
end

if length(sqzDims) == 1
    fids = zeros(dim_1, dim_2, dim_3, size(in.fids, dims.t));
    fids(1,1,1,:) = in.fids;
elseif length(sqzDims) == 2 
    if dims.extras==0 && dims.subSpecs==0 && dims.averages==0
        zeros(dim_1, dim_2, dim_3, size(in.fids, dims.t), size(in.fids, dims.coils));
        fids(1,1,1,:,:) = in.fids;
        dim_5_name = {'DIM_COIL'};
    elseif dims.extras==0 && dims.subSpecs==0 && dims.coils==0
        zeros(dim_1, dim_2, dim_3, size(in.fids, dims.t), size(in.fids, dims.averages));
        fids(1,1,1,:,:) = in.fids;
        dim_5_name = {'DIM_DYN'};
    elseif dims.extras==0 && dims.averages==0 && dims.coils==0
        zeros(dim_1, dim_2, dim_3, size(in.fids, dims.t), size(in.fids, dims.subSpecs));
        fids(1,1,1,:,:) = in.fids;
        dim_5_name = dim_names{1};
    end
   
end

% Initialize the NIfTI struct using the dicm2nii toolbox
% (https://github.com/xiangruili/dicm2nii)
try
    nii = nii_tool('init', fids);
catch ME
    switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            error(['Cannot find the function ''nii_tool.m''.' ...
                ' Please ensure that you have downloaded the required', ...
                ' dcm2nii toolbox (https://github.com/xiangruili/dicm2nii)', ...
                ' and added it to your MATLAB path.']);
        otherwise
            rethrow(ME);
    end
end

% Update NIfTI-MRS header
newDim = nii.hdr.dim;
newVoxOffset = nii.hdr.vox_offset;
nii.hdr = in.nii_mrs.hdr;
nii.img = conj(fids);

if exist('dim_5_name', 'var')
    nii.hdr_ext.dim_5 = dim_5_name;
    if exist('dim_6_name', 'var')
        nii.hdr_ext.dim_6 = dim_6_name;
        if exist('dim_7_name', 'var')
            nii.hdr_ext.dim_7 = dim_7_name;
        end
    end
end


nii.ext.ecode = 44;
nii.ext.edata_decoded = jsonencode(in.nii_mrs.hdr_ext);
len = int32(numel(nii.ext.edata_decoded));
myEdata = [nii.ext.edata_decoded'];
nii.ext.edata = myEdata;

% nii.hdr.pixdim     = in.nii_mrs.hdr.pixdim;
% nii.hdr.vox_offset     = in.nii_mrs.hdr.vox_offset;




% Save
nii_tool('save', nii, outfile);