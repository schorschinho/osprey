function osp_exportCSV (MRSCont,saveDestination)
%% osp_exportCSV (MRSCont,saveDestination)
%   This function creates CSV files from metabolite tables.
%
%   By default, all quantification tables which are initially created are
%   exported to the save destination QuantifyResults. The CSV files contain
%   3 columns. The first contains the concetration values, the second
%   contains a numerical code for the quantified metabolite, and the third
%   contains the name of the metabolite. This file can directly be used in
%   R (A R script will follow in the osprey/R folder soon)
%
%   USAGE:
%       MRSCont = OspreyQuantify(MRSCont,saveDestination);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       saveDestination = outputFolder.
%
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2019-11-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2019-11-06: First version of the code.

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