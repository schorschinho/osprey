function exportCSV (MRSCont,saveDestination,getResults)
if ~strcmp(MRSCont.opts.fit.method, 'LCModel')
    quants = {'amplMets','tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed'};
else
    quants = {'amplMets','tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed','CRLB','h2oarea'};
end
for ll = 1:length(getResults)
    for q = 1 : length(quants)
        if isfield(MRSCont.quantify.tables.(getResults{ll}), quants(q))
            if isstruct(MRSCont.quantify.tables.(getResults{ll}).(quants{q})) % To make export work on older MRSContainer
                if ~MRSCont.flags.isPRIAM
                    writetable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1,[saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_1.csv']);
                else
                    writetable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1,[saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_1.csv']);
                    writetable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_2,[saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_2.csv']);
                end
            else
                writetable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}),[saveDestination  filesep (getResults{ll}) '_' quants{q} '.csv']);
            end
        end
    end
end
end
