% fit_createMetabList.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% metabList = fit_createMetabList;
% 
% DESCRIPTION:
% Creates a list of metabolite basis functions that are to be included in 
% a fit.
% 
% OUTPUTS:
% metabList = structure including flags (1 = included, 0 = excluded) for
%             each metabolite included in the FID-A spin system definition,
%             plus various MM basis functions that may have been included
%             with fit_makeBasis.
%
% INPUTS:
% NONE

function metabList = fit_createMetabList

% Select metabolites to include in basis set
metabList.Ala      = 0;
metabList.Asc      = 0;
metabList.Asp      = 0;
metabList.bHB      = 0;
metabList.bHG      = 0;
metabList.Cit      = 0;
metabList.Cr       = 1;
metabList.CrCH2    = 1;
metabList.EtOH     = 0;
metabList.GABA     = 0;
metabList.GPC      = 0;
metabList.GSH      = 0;
metabList.Glc      = 0;
metabList.Gln      = 0;
metabList.Glu      = 0;
metabList.Gly      = 0;
metabList.H2O      = 1;
metabList.Ins      = 0;
metabList.Lac      = 0;
metabList.NAA      = 1;
metabList.NAAG     = 0;
metabList.PCh      = 0;
metabList.PCr      = 0;
metabList.PE       = 0;
metabList.Phenyl   = 0;
metabList.Scyllo   = 0;
metabList.Ser      = 0;
metabList.Tau      = 0;
metabList.Tyros    = 0;

% Select MM/lipid basis functions to include
metabList.MM09     = 1;
metabList.MM12     = 1;
metabList.MM14     = 1;
metabList.MM17     = 1;
metabList.MM20     = 1;
metabList.Lip09    = 1;
metabList.Lip13    = 1;
metabList.Lip20    = 1;

end