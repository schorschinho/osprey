% listValidBasisFunctionNames.m
% Georg Oeltzschner, Johns Hopkins University 2022
%
% USAGE:
% allMets = listValidBasisFunctionNames
% 
% DESCRIPTION:
% This function keeps a record of all valid metabolite names.
% If there is more than one common name for a metabolite, they are grouped
% inside another cell array. 
% The first name inside such a cell array is the default name for that
% metabolite in Osprey.
% 
% OUTPUTS:
% allMets = Cell array with a list of valid metabolite names.
%
% INPUTS
% type    = String (can be 'mets' or 'mm')

function listOfValidNames = listValidBasisFunctionNames(type)

switch type
    case 'mets'
        
        listOfValidNames =  {...
            'Ala',...   % Alanine
            'Asc',...   % Ascorbate
            'Asp',...   % Aspartate
            'bHB',...   % beta-hydroxybutyrate
            {'bHG', '2HG', '2-HG'},...   % 2-hydroxyglutarate
            'Cit',...   % Citrate
            {'Cr', 'Cre'},...    % Cr
            'Cystat',...% Cystat
            'CrCH2',... % negative CrCH2 correction signal
            'EA',...    % Ethanolamine
            'EtOH',...  % Ethanol
            'fCho',...  % free choline
            'GABA',...  % GABA
            'GPC',...   % Glycerophosphocholine
            'GSH',...   % Glutathione (reduced)
            'Glc',...   % Glucose
            'Gln',...   % Glutamine
            'Glu',...   % Glutamate
            {'Gly', 'Glyc'},...   % Glycine
            'H2O',...   % H2O
            {'mI', 'Ins', 'mIns'},...    % myo-inositol
            'Lac',...   % Lactate
            'NAA',...   % N-Acetylaspartate
            'NAAG',...  % N-Acetylaspartylglutamate
            'PCh',...   % Phosphocholine
            'PCr',...   % Phosphocreatine
            'PE',...    % Phosphoethanolamine
            'Phenyl',...    % Phenylalanine
            {'sI', 'Scyllo', 'sIns'},...    % scyllo-inositol
            'Ser',...   % Serine
            'Tau',...   % Taurine
            'Tyros',... % Tyrosine
            'NAA_Ace',...   % NAA acetyl
            'NAA_Asp',...   % NAA aspartyl
            };
        
    case 'mm'
        listOfValidNames = {...
            'MM09',...
            'MM12',...
            'MM14',...
            'MM17',...
            'MM20',...
            'MM22',...
            'MM27',...
            'MM30',...
            'MM32',...
            'Lip09',...
            'Lip13',...
            'Lip20',...
            'MM37',...
            'MM38',...
            'MM40',...
            'MM42',...
            {'MMexp','Mac'},... % Typical names for measured MM
            'MM_PRESS_PCC',...
            'MM_PRESS_CSO',...
            };

end

end


