function export_model_to_NII(MRSCont,data,ModelMatrix,outputFolder)
    
    non_zero = find(~cellfun('isempty', ModelMatrix));
    temp = ModelMatrix{non_zero(1)};
    metab_names = temp.BasisSets.names(logical(temp.BasisSets.includeInFit(end,:)));
    ModelMatrix=flip(ModelMatrix,1);
    if (MRSCont.raw{1}.nZvoxels > 1)
        ModelMatrix = flip(ModelMatrix,length(size(ModelMatrix)));
    end
    
    MRSI_model = get_MRSI_results(ModelMatrix);
    if (MRSCont.raw{1}.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
        shift = floor(data.nZvoxels/2);
        reorder = flip(1:data.nZvoxels);
        for ll = 1 : data.nZvoxels
                MaxData = max(max(max(real(squeeze(MRSI_model.data(:,:,:,ll))))));
                MinResidual = max(max(max(real(squeeze(MRSI_model.residual(:,:,:,ll))))));
    
                % MaxData = squeeze(max(real(squeeze(MRSI_model.data(:,:,:,ll))),[],1));
                % MinResidual = squeeze(min(real(squeeze(MRSI_model.residual(:,:,:,ll))),[],1));
                % MaxData = permute(repmat(MaxData,[1 1 size(MRSI_model.data,1)]),[3 1 2]);
                % MinResidual = permute(repmat(MinResidual,[1 1 size(MRSI_model.data,1)]),[3 1 2]);
                % 
                ToExport = data;
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.data(:,:,:,ll)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.data(:,:,:,ll));
                ToExport.nZvoxels = 1;
                ToExport.sz = size(ToExport.fids);
                ToExport.dims.averages = 0;
                ToExport.dims.subSpecs = 0;
                ToExport.dims.extras = 0;
                ToExport.dims.Xvoxels = 2;
                ToExport.dims.Yvoxels = 3;
                ToExport.dims.Zvoxels = 0;
                if isfield(ToExport.geometry,'slice_distance')
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                end
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit'))
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','data.nii.gz'));            
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.fit(:,:,:,ll)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.fit(:,:,:,ll));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','fit.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.baseline(:,:,:,ll)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.baseline(:,:,:,ll));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','baseline.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.residual(:,:,:,ll))+MaxData+abs(MinResidual),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.residual(:,:,:,ll))+MaxData+abs(MinResidual);
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','residual.nii.gz'));
    
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','metabs'))
                for mm = 1 : length(metab_names)
                    ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.metabs(:,mm,:,:,ll)),1), [], 1);
                    ToExport.specs = squeeze(MRSI_model.metabs(:,mm,:,:,ll));
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit','metabs',[metab_names{mm}  '.nii.gz']));
                end
    
            shift = shift - 1;
        end
    else
            % MaxData = max(real(squeeze(MRSI_model.data(:,:,:,:))),[],"all");
            % MinResidual = min(real(squeeze(MRSI_model.residual(:,:,:,:))),[],"all");

            MaxData = squeeze(max(real(squeeze(MRSI_model.data(temp.Data.ppm > temp.Options{temp.step}.optimFreqFitRange(1) & temp.Data.ppm < temp.Options{temp.step}.optimFreqFitRange(2),:,:,:,:))),[],1));
            MinResidual = squeeze(min(real(squeeze(MRSI_model.residual(temp.Data.ppm > temp.Options{temp.step}.optimFreqFitRange(1) & temp.Data.ppm < temp.Options{temp.step}.optimFreqFitRange(2),:,:,:,:))),[],1));
            
            if ndims(MaxData) == 4
                MaxData = permute(repmat(MaxData,[1 1 1 1 size(MRSI_model.data,1)]),[5 1 2 3 4]);
                MinResidual = permute(repmat(MinResidual,[1 1 1 1 size(MRSI_model.data,1)]),[5 1 2 3 4]);

                ToExport = data;
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.data(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.data(:,:,:,:,:));
                ToExport.sz = size(ToExport.fids);
    
    
                mkdir(fullfile(outputFolder,'fit'))
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','data.nii.gz'),{'DIM_EXP'});            
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.fit(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.fit(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','fit.nii.gz'),{'DIM_EXP'});
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.baseline(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.baseline(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','baseline.nii.gz'),{'DIM_EXP'});
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.residual(:,:,:,:,:))+(2*MaxData+abs(MinResidual)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.residual(:,:,:,:,:))+2*MaxData+abs(MinResidual);
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','residual.nii.gz'),{'DIM_EXP'});
    
                mkdir(fullfile(outputFolder,'fit','metabs'))
                for mm = 1 : length(metab_names)
                    ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.metabs(:,mm,:,:,:,:)),1), [], 1);
                    ToExport.specs = squeeze(MRSI_model.metabs(:,mm,:,:,:,:));
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','metabs',[metab_names{mm}  '.nii.gz']),{'DIM_EXP'});
                end
            end
            
            if ndims(MaxData) == 3
                MaxData = permute(repmat(MaxData,[1 1 1 size(MRSI_model.data,1)]),[4 1 2 3]);
                MinResidual = permute(repmat(MinResidual,[1 1 1 size(MRSI_model.data,1)]),[4 1 2 3]);
                ToExport = data;
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.data(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.data(:,:,:,:,:));
                ToExport.sz = size(ToExport.fids);
    
    
                mkdir(fullfile(outputFolder,'fit'))
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','data.nii.gz'));            
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.fit(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.fit(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','fit.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.baseline(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.baseline(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','baseline.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.residual(:,:,:,:,:))+(2*MaxData+abs(MinResidual)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.residual(:,:,:,:,:))+2*MaxData+abs(MinResidual);
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','residual.nii.gz'));
    
                mkdir(fullfile(outputFolder,'fit','metabs'))
                for mm = 1 : length(metab_names)
                    ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.metabs(:,mm,:,:,:,:)),1), [], 1);
                    ToExport.specs = squeeze(MRSI_model.metabs(:,mm,:,:,:,:));
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','metabs',[metab_names{mm}  '.nii.gz']));
                end
            end

            if ndims(MaxData) == 2
                MaxData = permute(repmat(MaxData,[1 1 size(MRSI_model.data,1)]),[3 1 2]);
                MinResidual = permute(repmat(MinResidual,[1 1 size(MRSI_model.data,1)]),[3 1 2]);
                ToExport = data;
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.data(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.data(:,:,:,:,:));
                ToExport.sz = size(ToExport.fids);
    
    
                mkdir(fullfile(outputFolder,'fit'))
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','data.nii.gz'));            
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.fit(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.fit(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','fit.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.baseline(:,:,:,:,:)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.baseline(:,:,:,:,:));
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','baseline.nii.gz'));
    
                ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.residual(:,:,:,:,:))+(2*MaxData+abs(MinResidual)),1), [], 1);
                ToExport.specs = squeeze(MRSI_model.residual(:,:,:,:,:))+2*MaxData+abs(MinResidual);
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','residual.nii.gz'));
    
                mkdir(fullfile(outputFolder,'fit','metabs'))
                for mm = 1 : length(metab_names)
                    ToExport.fids = ifft(ifftshift(squeeze(MRSI_model.metabs(:,mm,:,:,:,:)),1), [], 1);
                    ToExport.specs = squeeze(MRSI_model.metabs(:,mm,:,:,:,:));
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'fit','metabs',[metab_names{mm}  '.nii.gz']));
                end
            end

            


            
    end

end