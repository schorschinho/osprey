function exportQuant(MRSCont,saveDestination,getResults)

quants = {'amplMets','tCr','rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled','AlphaCorrWaterScaled','AlphaCorrWaterScaledGroupNormed'};
if strcmp(MRSCont.opts.fit.method, 'LCModel')
    quants = [quants, 'CRLB','h2oarea'];
end

for ll = 1:length(getResults)
    for q = 1 : length(quants)
        if isfield(MRSCont.quantify.tables.(getResults{ll}), quants(q))
            if isstruct(MRSCont.quantify.tables.(getResults{ll}).(quants{q})) % To make export work on older MRSContainer
                if ~MRSCont.flags.isPRIAM
                    MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1 = PopulateJSON(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, quants{q});
                    osp_WriteBIDsTable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, [saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_1']);
                else                    
                    MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1 = PopulateJSON(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, quants{q});
                    osp_WriteBIDsTable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, [saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_1']);
                    MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_2 = PopulateJSON(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, quants{q});
                    osp_WriteBIDsTable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_2, [saveDestination  filesep (getResults{ll}) '_' quants{q} '_Voxel_2']);
                end
            else
                MRSCont.quantify.tables.(getResults{ll}).(quants{q}) = PopulateJSON(MRSCont.quantify.tables.(getResults{ll}).(quants{q}).Voxel_1, quants{q});
                osp_WriteBIDsTable(MRSCont.quantify.tables.(getResults{ll}).(quants{q}), [saveDestination  filesep (getResults{ll}) '_' quants{q}])
            end
        end
    end
end


end

function[Table] = PopulateJSON(Table, Quant)

Table = addprop(Table, {'VariableLongNames'}, {'variable'}); % add long name to table properties

% Run thrrough quants to input relevant json info
switch Quant
    case'amplMets'
        L = 'Metabolite amplitude of %s';
        D = 'Raw amplitude of %s signal';
        u = 'arbitrary';
    case 'tCr'
        L = '%s to Cr ratio';
        D = 'Signal of %s represented as a ratio to the Cr signal';
        u = 'arbitrary';
    case 'rawWaterScaled'
        L = '%s to raw water ratio';
        D = 'Signal of %s represented as a ratio to the water peak';
        u = 'arbitrary';
    case 'CSFWaterScaled'
        L = '%s to raw water ratio, corrected for CSF fraction';
        D = 'Signal of %s as a ratio to raw water, but corrected for CSF fraction: WaterScaled/ (1 - fCSF)';
        u = 'arbitrary';
    case 'TissCorrWaterScaled'
        L = '%s molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006).';
        u = 'mol/Kg';
    case 'AlphaCorrWaterScaled'
        L = '%s alpha-corrected molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006). Alpha (cWM/cGM) correction applied according to: ConcIU/(fGM + alpha*fWM)';
        u = 'mol/Kg';
    case 'AlphaCorrWaterScaledGroupNormed'
        L = '%s alpha-corrected molal concentration, using tissue corrected water';
        D = '%s scaled to tissue corrected water according to Gasparovic et al, Magn Reson Med 55:1219-26 (2006). Alpha (cWM/cGM) correction applied according to: ConcIU/(fGM + alpha*fWM), with group normalization: (meanfGM + alpha * meanfWM) / ((fGM + alpha * fWM) * (meanfGM + meanfWM))';
        u = 'mol/Kg';
    case 'CRLB'
        L = '%s Cramer-Rao lower-bound';
        D = 'Lower bound of measurement error for %s';
        u = 'arbitrary';
    case 'h2oarea'
        L = 'Area ratio of %s to water peak';%CWDJ???
        D = 'Area of the water peak?';%CWDJ???
        u = 'arbitrary';
end

% input L, D, and u in to Table structure and parse this back to main function
for JJ=1:length(Table.Properties.VariableNames)
    Table.Properties.CustomProperties.VariableLongNames{JJ} = sprintf(L,Table.Properties.VariableNames{JJ});  
    Table.Properties.VariableDescriptions{JJ} = sprintf(D,Table.Properties.VariableNames{JJ});
    Table.Properties.VariableUnits{JJ} = u;
end

end