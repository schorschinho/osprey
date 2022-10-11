function osp_exportParams(MRSCont,path)
% osp_exportParam is used to export the fitting parameters optimized during 
% the Osprey fitting pipeline. The parameter fieldsnames are as follows:
%     - header      [metadata]
%     - ampl        
%     - ph0         [deg]
%     - ph1         [deg/ppm]
%     - gaussLB     [Hz]
%     - lorentzLB   [Hz]
%     - freqShift   [Hz]
%     - lineShape
%     - beta_j
%     - Res
%     - J
%     - refShift    [Hz]
%     - refFWHM     [Hz]
%     - ECC         [deg]
% 
% A header stores relevant metadata that describes the dataset corresponding 
% to the exportedparameters.
% 
  
% Define header fields
x = fieldnames(MRSCont.fit.resBasisSet.metab);
header.names         = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.name; 
header.nMets         = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.nMets; 
header.nMM           = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.nMM; 
header.spectralwidth = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.spectralwidth;
header.dwelltime     = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.dwelltime;
header.ppm           = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.ppm;
header.t             = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.t;
header.te            = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.te;
header.Bo            = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.Bo;
header.seq           = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.seq;
header.centerFreq    = MRSCont.fit.resBasisSet.metab.(x{1}){1,1}.centerFreq;
header.nDatasets     = MRSCont.nDatasets(1);


% Define the structure for storage
fields = fieldnames(MRSCont.fit.results.metab.fitParams{1, 1});
fields(strcmp(fields, 'prelimParams')) = [];
sz = [];
for i=1:length(fields)
    sz = [MRSCont.nDatasets(1), ...
          size(MRSCont.fit.results.metab.fitParams{1, 1}.(fields{i}))];
    params{i} = empty(sz);
end

% To add additional fitting parameters, the following two statements need to be 
% copied and modified to indicate the new parameter. Then add another line to 
% the end of the outer for-loop in the compiling section below.
params{length(params)+1} = empty([MRSCont.nDatasets(1), ...
                                  size(MRSCont.processed.metab{1,1}.phase_ecc)]);
fields{length(fields)+1} = 'ECC';


% Prepare for adding data
params = cell2struct(params, fields, dim=1);
params.header = header; % Add the header


% Compile the variables
for k = 1:MRSCont.nDatasets(1)
    for i = 1:length(fields)
        params.(fields{i})(k,:,:) = MRSCont.fit.results.metab.fitParams{1,k}(fields{i});
    end
    params.ECC(k,:,:) = MRSCont.processed.metab{1,k}.phase_ecc;
end


% Save the structure as a mat file
% The table option does not work because the contents do not have the same dimensions.
% writetable(array2table(params, 'VariableNames', names), 
%            fullfile(path,'parameters.csv'),'Delimiter',',');
save(fullfile(path,'parameters.mat'),'params')

  
end
