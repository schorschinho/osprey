function export_concs_to_NII(MRSCont,data,ModelMatrix,WaterModelMatrix,outputFolder, ModelMatrix2)
    if nargin == 6
        non_zero = find(~cellfun('isempty', ModelMatrix2));
        temp = ModelMatrix2{non_zero(1)};
        MRSI_model2 = get_MRSI_amplitudes(ModelMatrix2);
        MRSI_model2=flip(MRSI_model2,2);
        % if (MRSCont.raw_w{1}.nZvoxels > 1)
        %     MRSI_model2 = flip(MRSI_model2,length(size(MRSI_model2)));
        % end
    end

    non_zero = find(~cellfun('isempty', ModelMatrix));
    temp = ModelMatrix{non_zero(1)};
    metab_names = temp.Model{1}.CRLB.Properties.VariableNames;
    
    % if (MRSCont.raw{1}.nZvoxels > 1)
    %     ModelMatrix = flip(ModelMatrix,length(size(ModelMatrix)));
    % end

    shift = floor(data.nZvoxels/2);
    MRSI_model = get_MRSI_amplitudes(ModelMatrix);
    ModelMatrix=flip(ModelMatrix,2);
    if ~isempty(WaterModelMatrix)
        MRSI_model_water = get_MRSI_amplitudes(WaterModelMatrix);         
        MRSI_model_water.amplitudes = squeeze(sum(MRSI_model_water.amplitudes,1));
        MRSI_model_water.amplitudes = flip(MRSI_model_water.amplitudes ,1);
    end

    

    % if size(MRSI_model.amplitudes,1)==31 || size(MRSI_model.amplitudes,1)==28 %31 for all MMs
    %     metab_names(end+1:end+4) = {'tNAA','tCho','tCr','Glx'};
    % else
    %     metab_names(end+1:end+2) = {'tNAA','Glx'};
    % end

    if (MRSCont.raw{1}.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
        reorder = flip(1:data.nZvoxels);
        for ll = 1 : data.nZvoxels
                ToExport = data;
                if isfield(ToExport.geometry,'slice_distance')
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                end
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'fit'))
                out.hdr = ToExport.nii_mrs.hdr;
                out.hdr.dim(1) = 3;
                out.hdr.dim(2) = 1;
                out.hdr.pixdim(5) = 1;
                out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
                if ~isempty(WaterModelMatrix)
                    mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','rawWaterScaled'))
                end
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','amplitudes'))
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','tCr'))
                mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','CRLBs'))
                for mm = 1 : length(metab_names)
                    if ~isempty(WaterModelMatrix)
                        out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./MRSI_model_water.amplitudes(:,:,ll) *55500);         
                    
                        out.img(isnan(out.img)) =0;
                        out.img(isinf(out.img)) =0;
                        nii_tool('save', out, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','rawWaterScaled',[metab_names{mm}  '.nii.gz']));
                    end
                    if nargin == 5
                        out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./(squeeze(MRSI_model2.amplitudes(end-1,:,:,ll))));
                    else
                        out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./(squeeze(MRSI_model.amplitudes(end-1,:,:,ll))));
                    end
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','tCr',[metab_names{mm}  '.nii.gz']));
    
                    out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll)));
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','amplitudes',[metab_names{mm}  '.nii.gz']));
    
    
                    out.img = squeeze(MRSI_model.relCRLBs(mm,:,:,ll));
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','CRLBs',[metab_names{mm}  '_CRLBs.nii.gz']));
                end
                
    
            shift = shift - 1;
        end
    else
        ToExport = data;
        mkdir(fullfile(outputFolder,'fit'))
        out.hdr = ToExport.nii_mrs.hdr;
        out.hdr.dim(1) = 3;
        out.hdr.dim(2) = 1;
        out.hdr.pixdim(5) = 1;
        out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
        if ~isempty(WaterModelMatrix)
            mkdir(fullfile(outputFolder,'concs','rawWaterScaled'))
        end
        mkdir(fullfile(outputFolder,'concs','amplitudes'))
        mkdir(fullfile(outputFolder,'concs','tCr'))
        mkdir(fullfile(outputFolder,'concs','CRLBs'))
        for mm = 1 : length(metab_names)
            if ~isempty(WaterModelMatrix)
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./MRSI_model_water.amplitudes(:,:,:) *55500);         
            
                out.img(isnan(out.img)) =0;
                out.img(isinf(out.img)) =0;
                nii_tool('save', out, fullfile(outputFolder,'concs','rawWaterScaled',[metab_names{mm}  '.nii.gz']));
            end
            if nargin == 6
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./(squeeze(MRSI_model2.amplitudes(end-1,:,:,:))));
            else
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./(squeeze(MRSI_model.amplitudes(end-1,:,:,:))));
            end
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            nii_tool('save', out, fullfile(outputFolder,'concs','tCr',[metab_names{mm}  '.nii.gz']));

            out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:)));
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            nii_tool('save', out, fullfile(outputFolder,'concs','amplitudes',[metab_names{mm}  '.nii.gz']));


            out.img = squeeze(MRSI_model.relCRLBs(mm,:,:,:));
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            nii_tool('save', out, fullfile(outputFolder,'concs','CRLBs',[metab_names{mm}  '_CRLBs.nii.gz']));
        end
        out.img = squeeze(squeeze(MRSI_model_water.amplitudes(:,:,:)));
        out.img(isnan(out.img)) =0;
        out.img(isinf(out.img)) =0;
        nii_tool('save', out, fullfile(outputFolder,'concs','amplitudes',[  'water.nii.gz']));
    end

end
