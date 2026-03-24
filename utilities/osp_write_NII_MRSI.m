function nii = osp_write_NII_MRSI(in, outdir, UserNames, OspreyVersion, writeMask)

if nargin < 5
    writeMask = 0;
end
   
    % Get information about data 
    if in.nZvoxels > 1
        isMultiSlice = 1;
    else
        isMultiSlice = 0;
    end

    % Get temp geometry
    bckp_geometry = in.geometry;

    % Generate correct in-plane geometry from overall volume stored in SDAT/SPAR
    temp_geometry = in.geometry;
    temp_geometry.size.lr = temp_geometry.size.lr/in.nXvoxels;
    temp_geometry.size.ap = temp_geometry.size.ap/in.nYvoxels;

    % Generate correct slice thickness for multi slice MRSI data
    cc_org =  temp_geometry.size.cc;
    cc_size = cc_org - ((in.nZvoxels-1)*temp_geometry.slice_distance);
    cc_gap = temp_geometry.slice_distance - (cc_org-(in.nZvoxels-1)*temp_geometry.slice_distance);
    temp_geometry.slice_gap = cc_gap;
    temp_geometry.size.cc = cc_size;

    % Generate correct in-plane positions
    temp_geometry.pos.lr = temp_geometry.pos.lr + (in.nXvoxels/2 * temp_geometry.size.lr) - temp_geometry.size.lr;
    temp_geometry.pos.ap = temp_geometry.pos.ap + (in.nYvoxels/2 * temp_geometry.size.ap);

    % Generate correct slice location
    temp_geometry.pos.cc = temp_geometry.pos.cc - temp_geometry.size.cc/2 + cc_gap/2;

    if isMultiSlice
        if rem(in.nZvoxels, 2) == 1 %odd
            slice_shifts = (in.nZvoxels-1)/2 : -1 : -(in.nZvoxels-1)/2;
        else %even
            slice_shifts = (in.nZvoxels)/2 : -1 : -(in.nZvoxels)/2;
        end

        for zV = 1 : in.nZvoxels
            geometries{zV} =  temp_geometry;
            geometries{zV}.pos.cc = geometries{zV}.pos.cc + (slice_shifts(zV)*temp_geometry.slice_distance);
        end
    else
       geometries{1} =  temp_geometry; 
    end

    % Store data
    temp_fids = in.fids;
    temp_specs = in.specs;

    if ~writeMask
        for sl = 1 : length(geometries)
            in.geometry = geometries{sl};
            [~,filename,] = fileparts(in.OriginalFile);
            outfile = fullfile(outdir,[ filename '_vol_' num2str(sl)]);
            in.fids = squeeze(temp_fids(:,:,:,sl));
            in.specs = squeeze(temp_specs(:,:,:,sl));
            in.nZvoxels = 1;
            in.dims.Zvoxels =0;
            in.sz = size(in.fids);
            nii{sl} = io_writeniimrs(in, outfile, UserNames, OspreyVersion);
        end
    else
        for sl = 1 : length(geometries)
            in.geometry = geometries{sl};
            [~,filename,] = fileparts(in.OriginalFile);
            outfile = fullfile(outdir,[ filename '_vol_' num2str(sl)]);
            in.fids = abs(squeeze(temp_fids(:,:,:,sl)));
            in.specs = abs(squeeze(temp_specs(:,:,:,sl)));
            in.fids = in.fids(1,:,:);
            in.specs = in.specs(1,:,:);

            for y = 1 : in.nYvoxels
                for x = 1 : in.nXvoxels
                    index_mask_res(x,y) = str2num([sprintf('%02d',x) sprintf('%02d',y) ]);    
                end
            end
            in.fids(1,:,:) = index_mask_res;
            in.specs(1,:,:) = index_mask_res;

            in.nZvoxels = 1;
            in.dims.Zvoxels =0;
            in.sz = size(in.fids);
            nii{sl} = io_writeniimrs(in, outfile, UserNames, OspreyVersion);
        end
    end
end