function [out] = osp_SenseUnfolding(in,SENSE)

Voxels=2;
    %%

Navg = in.sz(in.dims.averages);

if in.dims.subSpecs ~= 0
    NsubSpec = in.sz(in.dims.subSpecs);
    fids = zeros(in.sz(in.dims.t),in.sz(in.dims.averages),in.sz(in.dims.subSpecs),Voxels);
    for ll = 1:Navg
        for mm = 1:NsubSpec
            fids(:,ll,mm,:) = (conj(SENSE.U) * squeeze(in.fids(:,:,ll,mm)).').';
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(fids(1,ll,mm,:)) ./ abs(conj(fids(1,ll,mm,:)));
            fids(:,ll,mm,:) = fids(:,ll,mm,:) .* permute(repmat(conj_norm(:,:).', [1 1 1 in.sz(in.dims.t) 1]), [4 2 3 1]);
        end
    end
else
    fids = zeros(in.sz(in.dims.t),in.sz(in.dims.averages),Voxels);
    for ll = 1:Navg
            fids(:,ll,:) = (conj(SENSE.U) * squeeze(in.fids(:,:,ll)).').';
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(fids(1,ll,:)) ./ abs(conj(fids(1,ll,:)));
            fids(:,ll,:) = fids(:,ll,:) .* permute(repmat(conj_norm(:,:).', [1 1 in.sz(in.dims.t) 1]), [3 2 1]);
    end
end

%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
if in.dims.t>in.dims.coils
    dims.t=in.dims.t-1;
else
    dims.t=in.dims.t;
end
dims.coils=0;
if in.dims.averages>in.dims.coils
    dims.averages=in.dims.averages-1;
else
    dims.averages=in.dims.averages;
end
if in.dims.subSpecs>in.dims.coils
    dims.subSpecs=in.dims.subSpecs-1;
else
    dims.subSpecs=in.dims.subSpecs;
end
if in.dims.extras>in.dims.coils
    dims.extras=in.dims.extras-1;
else
    dims.extras=in.dims.extras;
end

%re-calculate the sz variable
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;

%Adding MultiVoxelInfo
out.nXvoxels = Voxels;
out.nYvoxels = 0;
out.nZvoxels = 0;

%Adding voxel dimension
if in.dims.subSpecs ~= 0
    out.dims.Xvoxels = 4;
else
    out.dims.Xvoxels = 3;
end

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.MultiVoxel=1;
out.flags.addedrcvrs =1;


end
