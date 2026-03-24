function [MRSCont] = OspreyMRSIAtlasAnalysis(MRSCont)
    if ~isfield(MRSCont,'atlas')
        MRSCont.atlas = cell(1,1);
    end
    for kk = 1:MRSCont.nDatasets  
   
         [T1dir, T1name, T1ext]  = fileparts(MRSCont.files_nii{kk});
        if strcmp(T1ext,'.gz')
            T1name = strrep(T1name, '.nii','');
        end
        if strcmp(T1ext,'.gz')
            gunzip(MRSCont.files_nii{kk});
        end 


        switch MRSCont.opts.MRSI.atlas.name
            case 'AAL'
                atlas_path =  which('/libraries/AAL/modAAL3v1_1mm.nii'); 
                load(which('/libraries/AAL/modROI_MNI_V7_1mm_List.mat'));
                AtlasLabels = {ROI.Nom_L};
                ROIid = [ROI.ID];
                lr_index = cat(1,reshape([1:108],2,[])',reshape([117:138],2,[])',reshape([145:186],2,[])');
                tempAtlasLAbels = AtlasLabels(lr_index(:,1));
                tempAtlasLAbels = cellfun(@(s) s(1:end-2), tempAtlasLAbels, 'UniformOutput', false);
                tempAtlasLAbels{80} = [tempAtlasLAbels{80}(1:end-2) ' s'];
                tempAtlasLAbels{81} = [tempAtlasLAbels{80}(1:end-2) ' p'];
            case 'neuromorphometrics'
                atlas_path =  which('/libraries/neuromorphometrics/neuromorphometrics.nii');
                AtlasTable = readtable(which('/libraries/neuromorphometrics/neuromorphometrics.csv')); 
                ROIid = AtlasTable{:,1};
                AtlasLabels = AtlasTable{:,2};
                lr_index = cat(1,reshape([3:6],2,[])',reshape([8:15],2,[])',reshape([17:32],2,[])',reshape([37:136],2,[])');
                tempAtlasLAbels = AtlasLabels(lr_index(:,1));
                tempAtlasLAbels = cellfun(@(s) s(1:end-2), tempAtlasLAbels, 'UniformOutput', false);
        end
        warped_atlas_path = fullfile(T1dir,'Subject-Space-Atlas.nii');

        if ~isfield(MRSCont.atlas{kk},'fAtlas')
            % First we warp the atlas into subject space
            gunzip(fullfile(T1dir, ['iy_' T1name '.nii.gz'])); 
            warp_Atlas(fullfile(T1dir,['iy_' T1name '.nii']),atlas_path,fullfile(T1dir,[T1name '.nii']),warped_atlas_path);
        end

        Atlasvol  = spm_vol(warped_atlas_path);
        Atlasvol  = Atlasvol.private.dat(:,:,:);

        % Now we take the MRSI scan and calculate the lables for each MRSI voxel
        % based on a percentage threshold 

        vx =1; 
        total_voxels = MRSCont.raw{kk}.nXvoxels * MRSCont.raw{kk}.nYvoxels *MRSCont.raw{kk}.nZvoxels;
        index_mask = nii_tool('load', MRSCont.coreg.index_mask{kk}.fname);
        index_mask = double(index_mask.img);

        index_mask = round(index_mask);
        if ~isempty(MRSCont.coreg.gap_mask{kk})
            gap_mask = nii_tool('load', MRSCont.coreg.gap_mask{kk}.fname);
            gap_mask = double(gap_mask.img);
            gap_mask_logical = (gap_mask == 1);
        end

        % We need the number of labels 
        n_labels = length(ROIid);
        
        

        if ~isfield(MRSCont.atlas{kk},'fAtlas')

            MRSCont.atlas{kk}.fAtlas(:,:,:,:)  = zeros(n_labels,MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
    
            if isfield(MRSCont.opts.MRSI,'outerMask')
                if ~isfield(MRSCont.opts.MRSI.outerMask,'mask')  
                    MRSCont.opts.MRSI.outerMask.mask = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
                    MRSCont.opts.MRSI.outerMask.mask(MRSCont.opts.MRSI.outerMask.x(1):MRSCont.opts.MRSI.outerMask.x(2),...
                                                     MRSCont.opts.MRSI.outerMask.y(1):MRSCont.opts.MRSI.outerMask.y(2),...
                                                     MRSCont.opts.MRSI.outerMask.z(1):MRSCont.opts.MRSI.outerMask.z(2)) = 1;
                    % MRSCont.opts.MRSI.outerMask.mask = squeeze(MRSCont.opts.MRSI.outerMask.mask);
                end
            end
            for rr = 1 : MRSCont.raw{kk}.nZvoxels % This also loops over the different slice for MRSI
                for y = 1 : MRSCont.raw{kk}.nYvoxels
                    for x = 1 : MRSCont.raw{kk}.nXvoxels
                        if vx == 1    
                            msg = sprintf('Calculating anatomical label from voxel %3i out of %3i total voxels...\n', vx, total_voxels);
                            fprintf(msg);
                         else
                             msg = sprintf('Calculating anatomical label from voxel %3i out of %3i total voxels...\n', vx, total_voxels);
                                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                                fprintf([reverseStr, msg]);
                        end
                        if MRSCont.opts.MRSI.outerMask.mask(x,y,rr) && squeeze(MRSCont.seg.tissue.brain(kk,x,y,rr))
                            index = x * 1e6 + y * 1e3 + rr;
                            index_mask_temp =zeros(size(index_mask));
                            index_mask_temp(index_mask==index) =1;
    
                            % Apply gap mask if needed
                            if MRSCont.opts.MRSI.pseudo3D && ~isempty(MRSCont.coreg.gap_mask{kk})
                                index_mask_temp = index_mask_temp & ~gap_mask_logical;
                            end
    
                            mask_indices_temp = find(index_mask_temp);
                            
                            if ~isempty(mask_indices_temp)
                                % Extract atlas values only at mask locations
                                Atlas_vals = Atlasvol(mask_indices_temp);
                                n_voxels = length(mask_indices_temp);
                                
                                % Count occurrences of each ROI ID efficiently
                                for lb = 1:n_labels
                                    curr_id = ROIid(lb);
                                    
                                    % Count voxels with current ROI ID
                                    count = sum(Atlas_vals == curr_id);
                                    
                                    if count == 0
                                        fAtlas = 0;
                                    else
                                        fAtlas = count / n_voxels;
                                    end
                                    
                                    MRSCont.atlas{kk}.fAtlas(lb,x,y,rr) = fAtlas;
                                end
                            end
                        end
                        vx = vx + 1;                 
                    end
                end
            end
        end


        if ~exist(fullfile(MRSCont.outputFolder,['atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))]),'dir')
            mkdir(fullfile(MRSCont.outputFolder,['atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))]));
            % Export atlas in MRSI-space
            tempAtlas = MRSCont.atlas{kk}.fAtlas;
            tempAtlas(tempAtlas > MRSCont.opts.MRSI.atlas.AtlasThreshold) = 1;
            tempAtlas(tempAtlas < 1) = 0;
            ToExport = MRSCont.processed.A{kk};
            tempAtlas = permute(tempAtlas,[2 3 4 1]);
            out.hdr = ToExport.nii_mrs.hdr;
            out.hdr = ToExport.nii_mrs.hdr;
            out.hdr.dim(1) = 3;
            out.hdr.dim(2) = 1;
            out.hdr.dim(6) = n_labels;
            out.hdr.pixdim(5) = 1;
            out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes  
            out.img = tempAtlas;
            nii_tool('save', out, fullfile(MRSCont.outputFolder,['atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))],'MRSI-space_atlas.nii.gz'));  

            mkdir(fullfile(MRSCont.outputFolder,['atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))],'statistics'));
            mkdir(fullfile(MRSCont.outputFolder,['atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))],'pervoxel'));
        end

        

        % Run stat analysis 

        CombineLR = MRSCont.opts.MRSI.atlas.CombineLR;
        AtlasThreshold= MRSCont.opts.MRSI.atlas.AtlasThreshold;
        SNRThreshold = MRSCont.opts.MRSI.atlas.SNRThreshold;
        FWHMThreshold = MRSCont.opts.MRSI.atlas.FWHMThreshold;
        CRLBThreshold = MRSCont.opts.MRSI.atlas.CRLBThreshold;
        SDThreshold = MRSCont.opts.MRSI.atlas.SDThreshold;
        metabolites = MRSCont.opts.MRSI.atlas.metabolites;
        quantities = MRSCont.opts.MRSI.atlas.quantities;



        StatisticsAtlas = zeros(n_labels,8,length(metabolites),length(quantities));
        for qq = 1 : length(quantities)
            for lb = 1 : n_labels
                AtlasMask = squeeze(sum(MRSCont.atlas{kk}.fAtlas(lb,:,:,:),1));
                AtlasMaskIni = AtlasMask;

                 % Now lets apply other thresholds
                AtlasMask(MRSCont.quickMaps.A.FWHM > FWHMThreshold) = 0;
                AtlasMaskRemoveFWHM = AtlasMask;
                AtlasMask(MRSCont.quickMaps.A.SNR < SNRThreshold) = 0;
                AtlasMaskRemoveSNR = AtlasMask;

                % For CRLBs we need it per metabolite
                for mm = 1 : length(metabolites)
                    tempAtlasMask = AtlasMask;
                    tempAtlasMask(MRSCont.quantify.CRLBs.(metabolites{mm}) > CRLBThreshold(mm)) = 0;
                    AtlasMaskMetabs{mm} = tempAtlasMask;
                end

                AtlasMaskIni(AtlasMaskIni > AtlasThreshold) = 1;
                AtlasMaskIni(AtlasMaskIni < AtlasThreshold) = 0;
                AtlasMaskRemoveFWHM(AtlasMaskRemoveFWHM > AtlasThreshold) = 1;
                AtlasMaskRemoveFWHM(AtlasMaskRemoveFWHM < AtlasThreshold) = 0;
                AtlasMaskRemoveSNR(AtlasMaskRemoveSNR > AtlasThreshold) = 1;
                AtlasMaskRemoveSNR(AtlasMaskRemoveSNR < AtlasThreshold) = 0;

                for mm = 1 : length(metabolites)
                    tempAtlasMask = AtlasMaskMetabs{mm};
                    tempAtlasMask(tempAtlasMask > AtlasThreshold) = 1;
                    tempAtlasMask(tempAtlasMask < AtlasThreshold) = 0;
                    AtlasMaskMetabs{mm} = tempAtlasMask;
                end

                for mm = 1 : length(metabolites)
                    plotMap = MRSCont.quantify.(quantities{qq}).(metabolites{mm});

                    values = plotMap(AtlasMaskMetabs{mm} == 1);
                    mean_values = nanmean(values);
                    median_values = nanmedian(values);
                    std_values = nanstd(values);

                    ResultsAtlas{lb,mm,qq,:}=values;

                    if ~strcmp(quantities{qq},'CRLBs')
                        values(values > (mean_values + (SDThreshold*std_values))) = [];
                        values(values < (mean_values - (SDThreshold*std_values))) = [];
                        mean_values = nanmean(values);
                        median_values = nanmedian(values);
                        std_values = nanstd(values);
                    end

                    StatisticsAtlas(lb,1,mm,qq) = mean_values;
                    StatisticsAtlas(lb,2,mm,qq) = median_values;
                    StatisticsAtlas(lb,3,mm,qq) = std_values;
                    StatisticsAtlas(lb,4,mm,qq) = round(sum(AtlasMaskIni,'all'));
                    StatisticsAtlas(lb,5,mm,qq) = round(sum(AtlasMaskRemoveFWHM,'all'));
                    StatisticsAtlas(lb,6,mm,qq) = round(sum(AtlasMaskRemoveSNR,'all'));
                    StatisticsAtlas(lb,7,mm,qq) = round(sum(AtlasMaskMetabs{mm},'all'));
                    StatisticsAtlas(lb,8,mm,qq) = length(values);
                end

            end

            ind = n_labels + 1;


            % Let's quantify LR combined
            for lb = 1 : size(lr_index,1)
                AtlasMask = squeeze(sum(MRSCont.atlas{kk}.fAtlas(lr_index(lb,1):lr_index(lb,1),:,:,:),1));
                AtlasMaskIni = AtlasMask;

                 % Now lets apply other thresholds
                AtlasMask(MRSCont.quickMaps.A.FWHM > FWHMThreshold) = 0;
                AtlasMaskRemoveFWHM = AtlasMask;
                AtlasMask(MRSCont.quickMaps.A.SNR < SNRThreshold) = 0;
                AtlasMaskRemoveSNR = AtlasMask;

                % For CRLBs we need it per metabolite
                for mm = 1 : length(metabolites)
                    tempAtlasMask = AtlasMask;
                    tempAtlasMask(MRSCont.quantify.CRLBs.(metabolites{mm}) > CRLBThreshold(mm)) = 0;
                    AtlasMaskMetabs{mm} = tempAtlasMask;
                end

                AtlasMaskIni(AtlasMaskIni > AtlasThreshold) = 1;
                AtlasMaskIni(AtlasMaskIni < AtlasThreshold) = 0;
                AtlasMaskRemoveFWHM(AtlasMaskRemoveFWHM > AtlasThreshold) = 1;
                AtlasMaskRemoveFWHM(AtlasMaskRemoveFWHM < AtlasThreshold) = 0;
                AtlasMaskRemoveSNR(AtlasMaskRemoveSNR > AtlasThreshold) = 1;
                AtlasMaskRemoveSNR(AtlasMaskRemoveSNR < AtlasThreshold) = 0;

                for mm = 1 : length(metabolites)
                    tempAtlasMask = AtlasMaskMetabs{mm};
                    tempAtlasMask(tempAtlasMask > AtlasThreshold) = 1;
                    tempAtlasMask(tempAtlasMask < AtlasThreshold) = 0;
                    AtlasMaskMetabs{mm} = tempAtlasMask;
                end

                for mm = 1 : length(metabolites)
                    plotMap = MRSCont.quantify.(quantities{qq}).(metabolites{mm});

                    values = plotMap(AtlasMaskMetabs{mm} == 1);
                    mean_values = nanmean(values);
                    median_values = nanmedian(values);
                    std_values = nanstd(values);
                    ResultsAtlas{ind,mm,qq,:}=values;

                    if ~strcmp(quantities{qq},'CRLBs')
                        values(values > (mean_values + (SDThreshold*std_values))) = [];
                        values(values < (mean_values - (SDThreshold*std_values))) = [];
                        mean_values = nanmean(values);
                        median_values = nanmedian(values);
                        std_values = nanstd(values);
                    end

                    StatisticsAtlas(ind,1,mm,qq) = mean_values;
                    StatisticsAtlas(ind,2,mm,qq) = median_values;
                    StatisticsAtlas(ind,3,mm,qq) = std_values;
                    StatisticsAtlas(ind,4,mm,qq) = round(sum(AtlasMaskIni,'all'));
                    StatisticsAtlas(ind,5,mm,qq) = round(sum(AtlasMaskRemoveFWHM,'all'));
                    StatisticsAtlas(ind,6,mm,qq) = round(sum(AtlasMaskRemoveSNR,'all'));
                    StatisticsAtlas(ind,7,mm,qq) = round(sum(AtlasMaskMetabs{mm},'all'));
                    StatisticsAtlas(ind,8,mm,qq) = length(values);
                end
                ind = ind + 1;
            end    

        end
        
        AtlasLabels = [AtlasLabels,tempAtlasLAbels];      
        AtlasLabelsTable = table(AtlasLabels', 'VariableNames', {'Region'});
        MRSCont.AtlasResults{kk}.QuantitativeAALResults = StatisticsAtlas;
        MRSCont.AtlasResults{kk}.ResultsAAL = ResultsAtlas;
        % Now lets export
        ColNames = {'mean','median','std','total # voxels','# voxels FWHM filtered','# voxels FWHM + SNR filtered',...
                    '# voxels FWHM + SNR + CRLB filtered','# voxels FWHM + SNR + CRLB + SD filtered'};
        for qq = 1 : length(quantities)
            for mm = 1 : length(metabolites)
                % Write stats
                TableToExport = array2table(squeeze(StatisticsAtlas(:,:,mm,qq)),'VariableNames',ColNames);
                TableToExport = [AtlasLabelsTable, TableToExport];
                FileName = [MRSCont.outputFolder filesep 'atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))  filesep 'statistics' filesep 'AAL_statistics_' quantities{qq} '_' metabolites{mm}];
                writetable(TableToExport,[FileName '.txt'],'Delimiter','\t');
                movefile([FileName '.txt'],[FileName '.tsv']); % Change file extension to tsv

                % Write individual
                FileName = [MRSCont.outputFolder filesep 'atlas_' num2str(round(MRSCont.opts.MRSI.atlas.AtlasThreshold*100))  filesep 'pervoxel' filesep 'AAL_pervoxel_' quantities{qq} '_' metabolites{mm} '.tsv'];
                write_tsv_uneven(squeeze(ResultsAtlas(:,mm,qq))',FileName,AtlasLabels)
            end
        end
        delete(fullfile(T1dir, ['iy_' T1name '.nii'])); 
        gzip(warped_atlas_path);
        delete(warped_atlas_path);

    end
end

function warp_Atlas(y_field, aal_atlas, reference_img, output_file)
    % Warp with close parameters
    V_ref = spm_vol(reference_img);
    bb = spm_get_bbox(V_ref);
    vox = sqrt(sum(V_ref.mat(1:3,1:3).^2));
    % vox = abs(diag(V_ref.mat(1:3,1:3)))';
    
    % Add small buffer to ensure coverage
    bb_buffer = bb + [-1 -1 -1; 1 1 1];
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {y_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {aal_atlas};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb_buffer;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox * 0.99; % Slightly finer
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'temp_';
    
    spm_jobman('run', matlabbatch);
    
    % Force exact match using spm_reslice
    [pth, nm, ext] = fileparts(aal_atlas);
    temp_file = fullfile(pth, ['temp_' nm ext]);
    
    flags = struct('interp', 0, 'mask', 0, 'mean', 0, 'which', 1, 'prefix', '');
    spm_reslice({reference_img; temp_file}, flags);
    
    % Move file
    movefile(fullfile(pth, ['temp_' nm ext]), output_file);
end

function write_tsv_uneven(cellData,FileName,AtlasLabels)
    % Open file for writing
    fid = fopen(FileName, 'w');
    
    % Write header (optional)
    for lb = 1 : length(AtlasLabels)
        fprintf(fid, '%s\t', AtlasLabels{lb});
    end
    fprintf(fid, '\n');
    % Find the longest vector
    maxLen = max(cellfun(@length, cellData));
    
    % Write data row by row
    for row = 1:maxLen
    for col = 1:length(cellData)
        if row <= length(cellData{col})
            fprintf(fid, '%g', cellData{col}(row));
        end
        
        % Add tab if not the last column AND if there's more data to the right
        if col < length(cellData)
            % Check if current column has data OR any subsequent column has data for this row
            hasMoreData = false;
            for nextCol = col+1:length(cellData)
                if row <= length(cellData{nextCol})
                    hasMoreData = true;
                    break;
                end
            end
            
            if row <= length(cellData{col}) || hasMoreData
                fprintf(fid, '\t');
            end
        end
    end
    fprintf(fid, '\n');
end
    
    fclose(fid);
end