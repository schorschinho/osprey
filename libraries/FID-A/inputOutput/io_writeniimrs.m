% io_writeniimrs.m
% Georg Oeltzschner, Johns Hopkins University 2021.
% Adapted by
% Chris Davies-Jenkins, Johns Hopkins University 2022.
%
% USAGE:
% nii = io_writeniimrs(in, outfile, UserNames);
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

function nii = io_writeniimrs(in, outfile, UserNames, OspreyVersion)
%function nii = io_writeniimrs(in, outfile, UserNames)

% Bring FID array into NIfTI-MRS shape
dims = in.dims;

% Check InNames supplied for additional dimensions
if ~exist('UserNames','var')
    if dims.subSpecs || dims.extras
        error('Need InNames for additional dimensions')
    else
        UserNames = {};
    end
end

%% Manage essential dimensions
if ~isfield(dims, 'x') || dims.x==0
    dims.x = 0;
    dim_1  = 1;
else
    dim_1 = in.sz(dims.x);
end
if ~isfield(dims, 'y') || dims.y==0
    dims.y = 0;
    dim_2  = 1;
else
    dim_2 = in.sz(dims.y);
end
if ~isfield(dims, 'z')|| dims.z==0
    dims.z = 0;
    dim_3  = 1;
else
    dim_3 = in.sz(dims.z);
end

dim_4 = in.sz(dims.t);

ArraySize = [dim_1,dim_2,dim_3,dim_4];
dimorder = [dims.x,dims.y,dims.z,in.dims.t];
dimname = {};

%% Manage non-essential dimensions, starting with defaults
if dims.coils
    ArraySize = [ArraySize, in.sz(dims.coils)];
    dimorder = [dimorder,dims.coils];
    dimname = [dimname,{'DIM_COIL'}];
end
if dims.averages
    ArraySize = [ArraySize, in.sz(in.dims.averages)];
    dimorder = [dimorder,dims.averages];
    dimname = [dimname,{'DIM_DYN'}];
end

%% Now manage user defined dimensions
Cnt = 1;
if in.dims.subSpecs
    ArraySize = [ArraySize, in.sz(dims.subSpecs)];
    dimorder = [dimorder,dims.subSpecs];
    dimname = [dimname,UserNames(Cnt)];
    Cnt = Cnt+1;
end
if in.dims.extras
    ArraySize = [ArraySize, in.sz(dims.extras)];
    dimorder = [dimorder,dims.extras];
    dimname = [dimname,UserNames(Cnt)];
    Cnt = Cnt+1;
end

%% Do some error handling

if length(dimorder)>7
    error('%i dimensions detected, but maximum allowed is 7',length(dimorder))
elseif ~(length(UserNames)==(Cnt-1))
    error('%i user-defined dimensions, and %i dimension labels',Cnt-1,length(UserNames))
elseif length(UserNames)>3
    error('%i user-defined names supplied. Maximum of 3 allowed',length(UserNames))
end

%% Format fid output and export to nii

% If dimensions are ordered correctly, then just add singletons. Otherwise,
% reshape the array to correct for this.
%
% Note: Matlab will remove trailing singletons
if issorted(dimorder)
    fids = reshape(in.fids,ArraySize);
else
    %Permute fid non-singleton dimensions according to default BIDs precedence
    NonSingletonOrdered = permute(in.fids,nonzeros(dimorder));
    
    %Define overall array size including singleton dimensions
    ArraySize_S = double(ArraySize==1);
    ArraySize_S(~(ArraySize==1)) = size(NonSingletonOrdered);
    
    % Reshape output to include correct singletons 
    fids = reshape(NonSingletonOrdered,ArraySize_S);
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

%% Update NIfTI-MRS header
newDim = nii.hdr.dim;
newVoxOffset = nii.hdr.vox_offset;
% Add nii_mrs field if needed
if isfield(in,'nii_mrs')
    nii.hdr = in.nii_mrs.hdr;
else
    in.nii_mrs.hdr_ext = osp_generate_nii_hdr_ext(in,OspreyVersion);
    in.nii_mrs.hdr     = osp_generate_nii_hdr(in, nii.hdr,outfile);    
    nii.hdr = in.nii_mrs.hdr;
end
nii.img = conj(fids);

for JJ = 1:length(dimname)
    in.nii_mrs.hdr_ext.(sprintf('dim_%i',JJ+4)) = dimname{JJ};
end

if ~iscell(in.nii_mrs.hdr_ext.SpectrometerFrequency)
    in.nii_mrs.hdr_ext.SpectrometerFrequency = {in.nii_mrs.hdr_ext.SpectrometerFrequency};
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
end