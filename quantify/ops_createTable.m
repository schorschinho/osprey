function [MRSCont] = ops_createTable(MRSCont,qtfyType)

    % Extract metabolite names from basisset
    names = MRSCont.quantify.metabs;

    conc = zeros(MRSCont.nDatasets,length(names));
    for kk = 1:MRSCont.nDatasets
        conc(kk,:) = MRSCont.quantify.(qtfyType){kk}';
    end
    % Save back to Osprey data container
    MRSCont.quantify.tables.(qtfyType)  = array2table(conc,'VariableNames',names);
end

