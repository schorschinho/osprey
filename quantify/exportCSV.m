function exportCSV (MRSCont,saveDestination,getResults)
quants = {'tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed'};
for ll = 1:length(getResults)
    for q = 1 : length(fieldnames(MRSCont.quantify.tables.(getResults{ll})))
        if isfield(MRSCont.quantify.tables.(getResults{ll}), quants(q))
            writetable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}),[saveDestination  filesep (getResults{ll}) '_' quants{q} '.csv']);
        end
    end
end
end
