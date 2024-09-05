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
typicalMetaboliteCombinations = {'NAA','NAAG';'GPC','PCh';'Cr','PCr';'Glu','Gln';'EA','PE';'GABA','MM3co'; 'GABA', 'MM3to2'}; 
MetaboliteCombinationNames = {'tNAA','tCho','tCr','Glx','tEA','GABAplus','GABAplus'};


if ~strcmp(obj.Options{step}.parametrizations.metAmpl.type,'dynamic')
    nPars = 1;
    pos = 1;
else
    pos = find(contains(parametrizations.metAmpl.parameterNames,'Ampl'));   % For dynamic models we need to identify the reparameterized amplitude parameter position
    if isempty(pos)
        pos = 1;
    end
    nPars = length(parametrizations.metAmpl.parameterNames);
end
BMAT = zeros(size(invFisher,1),nPars*length(metaboliteNames));                          % Some LCModel nostalgia 
ll = 1;
for kk = parametrizations.metAmpl.start+(pos-1):nPars:parametrizations.metAmpl.end
   BMAT(kk,ll)=1;
   ll = ll + nPars;
end

%% 1. Get metabolite index, update Jacobian, update amplitude estimates
AddedMetaboliteCombinations = 0;
for mm = 1 : length(MetaboliteCombinationNames)
    idx_1 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,1}));        
    idx_2 = find(strcmp(metaboliteNames,typicalMetaboliteCombinations{mm,2}));
    if  ~isempty(idx_1) && ~isempty(idx_2)
        if ~strcmp(obj.Options{step}.parametrizations.metAmpl.type,'dynamic')
            AddedMetaboliteCombinations = AddedMetaboliteCombinations + 1;
            BMAT(parametrizations.metAmpl.start + idx_1 - 1,end+1)=1;
            BMAT(parametrizations.metAmpl.start + idx_2 - 1,end)=1;  
        else
            AddedMetaboliteCombinations = AddedMetaboliteCombinations + 1;
            BMAT(parametrizations.metAmpl.start + (nPars*idx_1 - nPars),end+1)=1;
            BMAT(parametrizations.metAmpl.start + (nPars*idx_2 - nPars),end)=1;  
        end
    end
end



%% 2. Calculate CRLB
if size(BMAT,2) > length(metaboliteNames)*nPars         % added new combinations 
    combinations = xk * BMAT;                           % build combinations of amplitude parameters
    DAPOSI = BMAT' * invFisher *BMAT;                   % multiply with inverse fisher matrix
    crlbs = sqrt(diag(DAPOSI));                         % get raw CRLBs values
    obj.Model{step}.Combined.rawCRLB.metAmpl = crlbs(end-AddedMetaboliteCombinations+1:end);  % Raw CRLBs for combined amplitudes                        
    relativeCRLB= (crlbs ./ combinations') * 100; % Relative CRLBs for combined amplitudes
    if nPars > 1
        relativeCRLB = cat(1,relativeCRLB(pos:nPars:end-AddedMetaboliteCombinations),relativeCRLB(end-AddedMetaboliteCombinations+1:end));
    end
    try
        obj.Model{step}.CRLB = array2table(relativeCRLB','VariableNames',[metaboliteNames MetaboliteCombinationNames(1:AddedMetaboliteCombinations)]'); % Save table with basis function names and relative CRLBs
    catch
    end  
end
end
