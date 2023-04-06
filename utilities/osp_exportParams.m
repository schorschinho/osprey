function osp_exportParams(MRSCont)
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
[hdr, label] = get_header(MRSCont);
header.names         = hdr.name; 
header.nMets         = hdr.nMets; 
header.nMM           = hdr.nMM; 
header.spectralwidth = hdr.spectralwidth;
header.dwelltime     = hdr.dwelltime;
header.ppm           = hdr.ppm;
header.t             = hdr.te;
header.Bo            = hdr.Bo;
header.seq           = hdr.seq;
header.centerFreq    = hdr.centerFreq;
header.nDatasets     = MRSCont.nDatasets(1);
header.isMRSI        = MRSCont.flags.isMRSI;


% Define the structure for storage
fitParams = get_fitParams(MRSCont);

num_voxels = numel(fitParams);
num_datasets = numel(fitParams{1}.(label{1}).fitParams); % ==MRSCont.nDatasets(1)
total = num_voxels * num_datasets;

fields = fieldnames(fitParams{1}.(label{1}).fitParams{1,1})';
% % % List fieldnames to ignore % % %
fields(strcmp(fields, 'prelimParams')) = [];
fields(strcmp(fields, 'LM_out')) = [];

for i=1:length(fields)
    sz = [total, size(fitParams{1}.(label{1}).fitParams{1}.(fields{i})), 1];
    params{i} = ones(sz);
end


% Used to create continuity between SVS and MRSI
MRSCont.processed = assertCell(MRSCont.processed);
MRSCont.QM        = assertCell(MRSCont.QM);
ref = strcmp(MRSCont.QM{1}.tables.Properties.VariableNames, "freqShift");
MRSCont.QM{1}.tables.Properties.VariableNames{ref} = char("global_freqShift");
qm_fields = MRSCont.QM{1}.tables.Properties.VariableNames; 


% To add additional fitting parameters, the following two statements need to be 
% copied and modified to indicate the new parameter. Then add another line to 
% the end of the outer for-loop in the compiling section below.
ignore = 0;
% % QM variables, i.e. Cr_SNR, Cr_FWHM, etc.
for i=1:length(qm_fields)
    params{length(params)+1} = ones([total, ...
                                     size(MRSCont.QM{1}.tables.(qm_fields{i})(1))]);
    fields{length(fields)+1} = qm_fields{i};
end

% % Eddy currents
if ismember('phase_ecc',fieldnames(MRSCont.processed{1}.(label{2}){1}))
    params{length(params)+1} = ones([total, ...
                                     size(MRSCont.processed{1}.(label{2}){1}.phase_ecc)]);
    fields{length(fields)+1} = 'ECC';
    ignore = ignore + 1;
    ecc = 1;
else
    ecc = 0;
end



% Prepare for adding data
params = cell2struct(params, fields, 2);
params.header = header; % Add the header

% Compile the variables
cnt = 0;
for m=1:num_voxels
    voxel = fitParams{m}.(label{1}).fitParams;
    for n=1:num_datasets
        cnt = cnt + 1;
        vxl = voxel{n};
        for i = 1:length(fields)-ignore
            if ~strcmp(fields{i},qm_fields)
                try
                    params.(fields{i})(cnt,:) = vxl.(fields{i});
                catch 
                    params.(fields{i})(cnt,:,:) = vxl.(fields{i});
                end
            else
                params.(fields{i})(cnt,:) = MRSCont.QM{1}.tables.(fields{i})(n,1);
            end
        end
        if ecc==1 
            params.ECC(cnt,:,:) = MRSCont.processed{m}.(label{2}){n}.phase_ecc;
        end
    end
end

% Save the structure as a mat file
% The table option does not work because the contents do not have the same dimensions.
% writetable(array2table(params, 'VariableNames', names), 
%            fullfile(path,'parameters.csv'),'Delimiter',',');
save(fullfile(MRSCont.opts.exportParams.path,
              'fitting_parameters.mat'),'-struct','params')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Auxiliary Functions    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hdr, label] = get_header(MRSCont)
    possible_keys = {{'metab','metab'},     % Osprey v2.0.0 and newer
                     {'off','A'},           % Osprey v1.2.0 and older (legacy)
                     {'diff1','diff1'},
                     {'diff2','diff2'}
                     {'sum','sum'}};          
    fields = fieldnames(MRSCont.fit.resBasisSet)';
    for i=1:length(fields)
        for n=1:length(possible_keys)
            if strcmpi(fields{i},possible_keys{n}{1})
                label = possible_keys{n};
            end
        end
    end
    if iscell(MRSCont.fit.resBasisSet.(label{1}))
        hdr = MRSCont.fit.resBasisSet.(label{1}){1,1};
    else % Assume struct
        x = fieldnames(MRSCont.fit.resBasisSet.(label{1}));
        hdr = MRSCont.fit.resBasisSet.(label{1}).(x{1}){1,1};
    end
end

function fitParams = get_fitParams(MRSCont)
    fitParams = MRSCont.fit.results;
    if isstruct(fitParams)
        fitParams = assertCell(fitParams);
    end
end

function cont = assertCell(cont)
    if ~iscell(cont)
        cont = {cont};
    end
end

end
