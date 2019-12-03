function exportCSV (MRSCont,saveDestination)
quants = {'tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed'};
for q = 1 : length(fieldnames(MRSCont.quantify.tables))
    if isfield(MRSCont.quantify.tables, quants(q))
        Concentration = table2array(MRSCont.quantify.tables.(quants{q}));
        names = cell(size(Concentration));
        Concentration = reshape(Concentration,[],1);
        Concentration(:,2) = ones(MRSCont.nDatasets*length(MRSCont.quantify.metabs),1);    
        for m = 0 : length(MRSCont.quantify.metabs)-1
        Concentration((1+MRSCont.nDatasets*m):(MRSCont.nDatasets*(m+1)),2) = ones(MRSCont.nDatasets,1)* (m+1);    
            for kk = 1 : MRSCont.nDatasets
                names{kk,m+1} = MRSCont.quantify.metabs{m+1};
            end        
        end
        names = reshape(names,[],1);
        ConC = Concentration(:,1);
        MetaboliteNumber = Concentration(:,2);
        MetaboliteString = names;

        data_table = table(ConC,MetaboliteNumber,MetaboliteString);
        writetable(data_table,[saveDestination  filesep quants{q} '.csv']);
    end
end
end