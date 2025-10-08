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
            'AcAc',...  % Acetoacetate
            {'Ace', 'Act'},...   % Acetate
            {'AcO', 'Acn'},...   % Acetone
            'Ala',...   % Alanine
            'Asc',...   % Ascorbate
            'Asp',...   % Aspartate
            'Bet',...   % Betaine
            {'bHB', 'bHb'},...   % beta-hydroxybutyrate
            {'bHG', '2HG', '2-HG'},...   % 2-hydroxyglutarate
            'Car', ...  % Carnitine
            'Cit',...   % Citrate
            {'Cr', 'Cre'},... % Cr
            'Cr_methyl_only',...% Cr_methyl_only
            'Cr_methylene_only',...%Cr_methylene_only
            'Cys', ...  % Cysteic acid
            'Cystat',...% Cystat
            'CrCH2',... % negative CrCH2 correction signal
            'EA',...    % Ethanolamine
            {'EtOH', 'Eth'},...  % Ethanol
            {'fCho', 'Cho'},...  % free choline
            'Fuc',...   % Fucose
            'GABA',...  % GABA
            'Gcn',...   % Glucone
            'Gcr',...   % Glucoronic acid
            'GPC',...   % Glycerophosphocholine
            'GSH',...   % Glutathione (reduced)
            'Glc',...   % Glucose
            'Gln',...   % Glutamine
            'Glu',...   % Glutamate
            {'Gly', 'Glyc'},...   % Glycine
            'Gua',...   % Guanidinoacetate
            'H2O',...   % H2O
            'HCar',...  % Homocarnosine
            {'ILc','iLe','Ile'},...   % Isoleucine
            {'mI', 'Ins', 'mIns'},...    % myo-inositol
            'Lac',...   % Lactate
            'Leu',...   % Leucine
            'Lys',...   % Lysine
            'MSM',...   % Methyl-sulfonyl-methane
            'NAA',...   % N-Acetylaspartate
            'NAA_Acetyl_only',...
            'NAA_Aspartate_only',...   
            'NAAG',... % N-Acetylaspartylglutamate
            'NAAG_Acetyl_only',...
            'NAAG_AspGlu_only',... 
            {'PCh', 'PCho'},...   % Phosphocholine
            'PCr',...   % Phosphocreatine
            'PE',...    % Phosphoethanolamine
            'Pgc',...   % Propyleneglycol
            {'Phenyl', 'PAl','Phe'},...    % Phenylalanine
            'Pyr',...   % Pyruvate
            {'sI', 'Scyllo', 'sIns'},...    % scyllo-inositol
            'Ser',...   % Serine
            'Suc',...   % Succinate
            'Tau',...   % Taurine
            'Thr',...   % Threonine
            {'Tryp','Trp'},...  % Tryptophan
            {'Tyros','Tyr'},... % Tyrosine
            'Val',...   % Valine
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
            {'MMexp','Mac', 'MMmeas'},... % Typical names for measured MM
            'MM_PRESS_PCC',...
            'MM_PRESS_CSO',...
            };

end

end


