function obj = calculateCombinedCRLB(obj, invFisher, xk, metaboliteNames, parametrizations)
%%  [obj] = calculateCombinedCRLB(jacobian, metaboliteNames, parametrizations)
%   This method calculates the CRLBs of metabolite combinations inthe OspreyFitObj 
%
%   USAGE:
%       [obj] = calculateCombinedCRLB(obj,jac, metaboliteNames, parametrizations)
%
%   INPUTS:
%       invFisher   = inverse fisher matrix
%       xk          = final parameter vector
%       metaboliteNames  = metabolites in model
%       parametrizations  = struct with parametrizations in xk vector
%       
%   OUTPUTS:
%       obj       = OspreyFitObj with updated parametrizations.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)

%% 0. Set names cell arrays and houskeeping
step            = obj.step;                                                                     % Get step counter
typicalMetaboliteCombinations = {'NAA','NAAG';'GPC','PCh';'Cr','PCr';'Glu','Gln';'EA','PE';...
                                 'NAA_Acetyl_only','NAAG_Acetyl_only';'GPC_pCh2_only','PCh_trimethyl_only';
                                 'Cr_methyl_only','PCr_ch3_only';'Cr_methylene_only','PCr_ch2nhnh_only'}; 
MetaboliteCombinationNames = {'tNAA','tCho','tCr','Glx','tEA','tNAA_Acetyl_only','tCho_pCh2_only','tCr_methyl_only','tCr_mehtylene_only'};

BMAT = zeros(size(invFisher,1),length(metaboliteNames));                          % Some LCModel nostalgia 
ll = 1;
for kk = parametrizations.metAmpl.start:parametrizations.metAmpl.end
   BMAT(kk,ll)=1;
   ll = ll + 1;
end
%% 1. Get metabolite index, update Jacobian, update amplitude estimates
AddedMetaboliteIndex =[];
AddedMetaboliteCombinations = 0;
for mm = 1 : length(MetaboliteCombinationNames)
    idx_1 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,1}));        
    idx_2 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,2}));
    if  ~isempty(idx_1) && ~isempty(idx_2)
        AddedMetaboliteCombinations = AddedMetaboliteCombinations + 1;
        BMAT(parametrizations.metAmpl.start + idx_1 - 1,end+1)=1;
        BMAT(parametrizations.metAmpl.start + idx_2 - 1,end)=1;    
        AddedMetaboliteIndex(end+1) = mm;
    end
end



%% 2. Calculate CRLB
if size(BMAT,2) > length(metaboliteNames)               % added new combinations 
    combinations = xk * BMAT;                           % build combinations of amplitude parameters
    DAPOSI = BMAT' * invFisher *BMAT;                   % multiply with inverse fisher matrix
    crlbs = sqrt(diag(DAPOSI));                         % get raw CRLBs values
    obj.Model{step}.Combined.rawCRLB.metAmpl = crlbs(end-AddedMetaboliteCombinations+1:end);  % Raw CRLBs for combined amplitudes                        
    relativeCRLB= (crlbs ./ combinations') * 100; % Relative CRLBs for combined amplitudes
    
    try
        obj.Model{step}.CRLB = array2table(relativeCRLB','VariableNames',[metaboliteNames MetaboliteCombinationNames(AddedMetaboliteIndex)]'); % Save table with basis function names and relative CRLBs
    catch
    end  
end
end
