function MRSI_model = get_MRSI_amplitudes(FitMatrix)
    dims = size(FitMatrix);
    if length(dims) == 2
        dims(end+1) = 1;
    end
    non_zero = find(~cellfun('isempty', FitMatrix));
    temp = FitMatrix{non_zero(1)};

    metaboliteNames = temp.Model{end}.CRLB.Properties.VariableNames;

    MRSI_model.amplitudes = zeros(size(temp.Model{end}.CRLB,2),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3));
    MRSI_model.relCRLBs = zeros(size(temp.Model{end}.CRLB,2),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3));

    typicalMetaboliteCombinations = {'NAA','NAAG';'GPC','PCh';'Cr','PCr';'Glu','Gln';'EA','PE';...
                                 'NAA_Acetyl_only','NAAG_Acetyl_only';'GPC_pCh2_only','PCh_trimethyl_only';
                                 'Cr_methyl_only','PCr_ch3_only';'Cr_methylene_only','PCr_ch2nhnh_only'}; 
    MetaboliteCombinationNames = {'tNAA','tCho','tCr','Glx','tEA','tNAA_Acety_only','tCho_pCh2_only','tCr_methyl_only','tCr_mehtylene_only'};

    AddedMetaboliteIndex =[];
    AddedMetaboliteIndices = [];
    AddedMetaboliteCombinations = 0;
    for mm = 1 : length(MetaboliteCombinationNames)
        idx_1 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,1}));        
        idx_2 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,2}));
        if  ~isempty(idx_1) && ~isempty(idx_2)                                   
            AddedMetaboliteCombinations = AddedMetaboliteCombinations + 1;
            AddedMetaboliteIndex(end+1) = mm;
            AddedMetaboliteIndices(1,end+1) = idx_1;
            AddedMetaboliteIndices(2,end) = idx_2;
        end
    end

    for z = 1: dims(3)
        for x = 1: dims(1)
            for y = 1: dims(2)
                if ~isempty(FitMatrix{x,y,z})
                    MRSI_model.amplitudes(1:size(temp.Model{end}.fit.metabs,2),x,y,z) = FitMatrix{x,y,z}.Model{end}.parsOut.metAmpl(1,:);   
                    for mm = 1 : AddedMetaboliteCombinations
                        MRSI_model.amplitudes(size(temp.Model{end}.fit.metabs,2)+mm,x,y,z) = FitMatrix{x,y,z}.Model{end}.parsOut.metAmpl(1,AddedMetaboliteIndices(1,mm)) + FitMatrix{x,y,z}.Model{end}.parsOut.metAmpl(1,AddedMetaboliteIndices(2,mm));
                    end                    
                    MRSI_model.relCRLBs(:,x,y,z) = FitMatrix{x,y,z}.Model{end}.CRLB{:,:};
                end
            end
        end
    end

end