function [vol_mask, T1_max, voxel_ctr,index_mask,gap_mask] = coreg_MRSI(in,vol_image,outputFolder,pseudo3D,nii_shifts)
    
    outputFolder = fullfile(outputFolder,'VoxelMasks');
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    gap_mask = [];
    %% This is voxel masks with ones
    % We can use it for a grid
   
    if (in.nZvoxels > 1) && ~pseudo3D 
        shift = floor(in.nZvoxels/2);
        reorder = flip(1:in.nZvoxels);
        for ll = 1 : in.nZvoxels
                ToExport = in;
                if isfield(ToExport.geometry,'slice_distance')
                    VoxelShift = [-in.nXvoxels/2 + nii_shifts(1) , -in.nYvoxels/2 + nii_shifts(2), ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                end
                out.hdr = ToExport.nii_mrs.hdr;
                [out] = updateNiiHDR(out);

                out.img = ones(size(squeeze(in.fids(1,:,:,ll))));
                nii_tool('save', out, fullfile(outputFolder,['VoxelMask_slice_' num2str(reorder(ll)) '.nii.gz']));               
                shift = shift - 1;
        end
    else
        VoxelMask.hdr = in.nii_mrs.hdr;
        [VoxelMask] = updateNiiHDR(VoxelMask);  
        if in.dims.subSpecs == 0
            VoxelMask.img = ones(size(squeeze(in.fids(1,:,:,:,:))));
        else
            if in.nZvoxels > 1
                VoxelMask.img = ones([in.nXvoxels,in.nYvoxels,in.nZvoxels]);
            else
                VoxelMask.img = ones([in.nXvoxels,in.nYvoxels]);
            end
        end
        nii_tool('save', VoxelMask, fullfile(outputFolder,'VoxelMask.nii.gz'));   
        
        % We want the gap slices as separate nifti files here
        if isfield(in.geometry,'gap')
            shift = floor(in.nZvoxels/2);
            reorder = flip(1:in.nZvoxels);
            shift = 1;
            for ll = 1 : in.nZvoxels
                    ToExport = in;                   

                    ToExport.geometry.size.cc = ToExport.geometry.gap/2;
                    ToExport.nii_mrs.hdr     = osp_generate_nii_hdr(ToExport, ToExport.nii_mrs.hdr,ToExport.OriginalFile); 
                    VoxelShift = [-in.nXvoxels/2 + nii_shifts(1) , -in.nYvoxels/2 + nii_shifts(2), (in.geometry.slice_distance/ToExport.geometry.size.cc)*shift + floor((in.geometry.size.cc-in.geometry.gap)/ToExport.geometry.size.cc/2) + 1];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice               
                    out.hdr = ToExport.nii_mrs.hdr;
                    [out] = updateNiiHDR(out);
                    if ToExport.dims.subSpecs == 0
                        out.img = ones(size(squeeze(in.fids(1,ToExport.nXvoxels,ToExport.nYvoxels,ll))));
                    else
                        out.img = ones(1,ToExport.nXvoxels,ToExport.nYvoxels);
                    end
                    nii_tool('save', out, fullfile(outputFolder,['VoxelMask_slice_' num2str(reorder(ll)) '_gap_1.nii']));  

                    ToExport = in;
                    ToExport.geometry.size.cc = ToExport.geometry.gap/2;
                    ToExport.nii_mrs.hdr     = osp_generate_nii_hdr(ToExport, ToExport.nii_mrs.hdr,ToExport.OriginalFile);                  
                    VoxelShift = [-in.nXvoxels/2 + nii_shifts(1) , -in.nYvoxels/2 + nii_shifts(2), -1*((in.geometry.slice_distance/ToExport.geometry.size.cc)*shift + floor((in.geometry.size.cc-in.geometry.gap)/ToExport.geometry.size.cc/2) + 1)];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                    out.hdr = ToExport.nii_mrs.hdr;
                    [out] = updateNiiHDR(out);
    
                    if ToExport.dims.subSpecs == 0
                        out.img = ones(size(squeeze(in.fids(1,ToExport.nXvoxels,ToExport.nYvoxels,ll))));
                    else
                        out.img = ones(1,ToExport.nXvoxels,ToExport.nYvoxels);
                    end
                    nii_tool('save', out, fullfile(outputFolder,['VoxelMask_slice_' num2str(reorder(ll)) '_gap_2.nii']));
                    shift = shift - 1;
            end
        end
    end

    
    % Now let's generate a grid
    gunzip(fullfile(outputFolder,'VoxelMask.nii.gz'));
    nii_file = fullfile(outputFolder,'VoxelMask.nii');
    vol_mask = spm_vol(fullfile(outputFolder,'VoxelMask.nii'));

    V = spm_vol(nii_file);
    Y = spm_read_vols(V);

    % V.mat = [VoxelMask.hdr.srow_x; VoxelMask.hdr.srow_y; VoxelMask.hdr.srow_z;];
    % Y = VoxelMask.img;
    
    % Get voxel size and origin
    voxel_size = sqrt(sum(V.mat(1:3,1:3).^2));
    origin = V.mat(1:3,4);

    rotation_scaling_matrix = V.mat(1:3, 1:3);
    normalized_rotation = rotation_scaling_matrix ./ voxel_size;  % Normalize rotation
    
    % Find non-zero voxels
    [idx_x, idx_y, idx_z] = ind2sub(size(Y), find(Y > 0));
    num_voxels = length(idx_x);
    
    % Preallocate mesh data
    all_vertices = [];
    all_faces = [];
    face_template = [
        1 2 4; 1 4 3;  % Bottom
        5 6 8; 5 8 7;  % Top
        1 2 6; 1 6 5;  % Side 1
        2 4 8; 2 8 6;  % Side 2
        4 3 7; 4 7 8;  % Side 3
        3 1 5; 3 5 7   % Side 4
    ];
    v0 = [0 0 0; 1 0 0; 0 1 0; 1 1 0;
          0 0 1; 1 0 1; 0 1 1; 1 1 1];
    
    v0 = (v0 - 0.5);  % center at origin
    
    % Generate cubes for each voxel
    for i = 1:num_voxels
        % Voxel center in mm
        ijk = [idx_x(i), idx_y(i), idx_z(i), 1]';
        xyz = V.mat * ijk;
        
        % Scale and shift cube vertices
        verts = (v0 .* voxel_size * normalized_rotation')  + xyz(1:3)';
        % verts = (v0 .* voxel_size) * V.mat(1:3,1:3)' + xyz(1:3)';
        
        % Add to global list
        offset = size(all_vertices, 1);
        all_vertices = [all_vertices; verts];
        all_faces = [all_faces; face_template + offset];
    end
    
    % Create GIfTI surface structure
    g = gifti;
    g.faces = int32(all_faces);
    g.vertices = single(all_vertices);
    
    gii_filename_VoxelGrid = fullfile(outputFolder,'VoxelMaskGrid.gii');
    save(g, gii_filename_VoxelGrid);

    delete(fullfile(outputFolder,'VoxelMask.nii'))
    %% Index mask 
    % This has indices for all voxel positions
    if in.dims.subSpecs == 0
        index_mask = zeros(size(squeeze(in.fids(1,:,:,:))));
    else
        if in.nZvoxels > 1
            index_mask = zeros([in.nXvoxels,in.nYvoxels,in.nZvoxels]);
        else
           index_mask = zeros([in.nXvoxels,in.nYvoxels]);
        end
    end
    for z = 1 : size(index_mask,3)
        for y = 1 : size(index_mask,2)
            for x = 1 : size(index_mask,1)
                    index_mask(x,y,z) = str2num([sprintf('%03d',x) sprintf('%03d',y) sprintf('%03d',z) ]);    
            end
        end
    end

    if (in.nZvoxels > 1) && ~pseudo3D 
        reorder = flip(1:data.nZvoxels);
        for ll = 1 : data.nZvoxels
                ToExport = in;
                if isfield(ToExport.geometry,'slice_distance')
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                end
                out.hdr = ToExport.nii_mrs.hdr;
                [out] = updateNiiHDR(out);

                out.img = squeeze(index_mask(:,:,ll));
                nii_tool('save', out, fullfile(outputFolder,['IndexMask_slice_' num2str(reorder(ll)) '.nii.gz']));
                gunzip(fullfile(outputFolder,['IndexMask_slice_' num2str(reorder(ll)) '.nii.gz']));
                shift = shift - 1;
        end
    else
        ToExport = in;
        out.hdr = ToExport.nii_mrs.hdr;
        out.hdr.dim(1) = 3;
        out.hdr.pixdim(1) = -1;
        out.img = index_mask;
        nii_tool('save', out, fullfile(outputFolder,'IndexMask.nii.gz'));   
        gunzip(fullfile(outputFolder,'IndexMask.nii.gz'));
    end

    %% Reslice the index masks and gap masks to T1 image

    % Create SPM volume and read in the NIfTI file with the structural image.
    [T1,XYZ]    = spm_read_vols(vol_image);
    T1_max      = max(T1(:));

    geom = in.geometry;
    [BB,vx] = spm_get_bbox(vol_image.fname);
    voxel_ctr(:,:,:) = [0 0 0];
    resize_img(fullfile(outputFolder,'IndexMask.nii'), vol_image, 1);
    index_mask = spm_vol(fullfile(outputFolder, 'rIndexMask.nii'));
    gzip(fullfile(outputFolder, 'rIndexMask.nii'));
    delete(fullfile(outputFolder,'IndexMask.nii'))
    delete(fullfile(outputFolder,'rIndexMask.nii'))

    if isfield(in.geometry,'gap')  && pseudo3D 
        reslice_gap_files = {};
        for ll = 1 : in.nZvoxels
            resize_img(fullfile(outputFolder,['VoxelMask_slice_' num2str(ll) '_gap_1.nii']), vol_image, 1);
            delete(fullfile(outputFolder,['VoxelMask_slice_' num2str(ll) '_gap_1.nii']));
            reslice_gap_files{end+1} = fullfile(outputFolder,['rVoxelMask_slice_' num2str(ll) '_gap_1.nii']);
            resize_img(fullfile(outputFolder,['VoxelMask_slice_' num2str(ll) '_gap_2.nii']), vol_image, 1);
            delete(fullfile(outputFolder,['VoxelMask_slice_' num2str(ll) '_gap_2.nii']));
            reslice_gap_files{end+1} = fullfile(outputFolder,['rVoxelMask_slice_' num2str(ll) '_gap_2.nii']);
        end
        vol_gap = spm_vol(reslice_gap_files{1});
        % for gg = 2 : length(reslice_gap_files)
        %     gap_vol_temp = spm_vol(reslice_gap_files{gg});
        %     vol_gap.private.dat(:,:,:) = vol_gap.private.dat(:,:,:) + gap_vol_temp.private.dat(:,:,:);
        % end
        n_files = length(reslice_gap_files);
        all_data = zeros([vol_gap.dim, n_files]);
        for gg = 1:n_files
            gap_vol_temp = spm_vol(reslice_gap_files{gg});
            all_data(:,:,:,gg) = gap_vol_temp.private.dat(:,:,:);
            delete(reslice_gap_files{gg});
        end
        vol_gap.private.dat(:,:,:) = sum(all_data, 4);
        gap_vol      = vol_gap.private.dat(:,:,:);
        vol_gap.fname    = fullfile(outputFolder, 'GapMask.nii');
        vol_gap.descrip  = ['MRSI_gap'];
        vol_gap          = spm_write_vol(vol_gap, gap_vol);

        gap_mask = spm_vol(fullfile(outputFolder, 'GapMask.nii'));
        delete(reslice_gap_files{1});
        gzip(fullfile(outputFolder, 'GapMask.nii'));
        delete(fullfile(outputFolder,'GapMask.nii'));
    end
    

    % realignT1MRSI(vol_image.fname,fullfile(outputFolder, 'rIndexMask.nii'),gap_mask.fname)
end



function realignT1MRSI(T1file,index_mask,gap_mask)

matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                                                    {[T1file ',1']},
                                                    {[index_mask ',1']},
                                                    {[gap_mask ',1']},
                                                    }';
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'a';

spm_jobman('run',matlabbatch);

end

function [out] = updateNiiHDR(in)
    out = in;
    out.hdr.dim(1) = 3;
    out.hdr.dim(2) = 1;
    out.hdr.pixdim(5) = 1;
    out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes  
    if out.hdr.dim(3) > 1
        out.hdr.dim(3) = 1;
    end
end