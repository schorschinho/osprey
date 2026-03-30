function [MRSCont] = OspreyMRSIOverview(MRSCont)
%% Run additional analysis with Osprey Overview MRSI
%% Representative Spectra
% Calculate representative white matter and gray matter spectra according
% to Gorywala et al. MRM 2017.
%
% Input options include slices to pick and the WM + GM threshhold to pick
% spectra for the calculation.

MRSCont = computeRepresentativeSpectra(MRSCont,'A');
%% Global Concentration with linear regression

%   This function calculates global concentrations using a linear
%   regression appoach as drescribed in Tal A, Kirov II, Grossman RI, Gonen O. 
%   The role of gray and white matter segmentation in quantitative proton MR 
%   spectroscopic imaging. NMR Biomed. 2012;25(12):1392-1400. doi:10.1002/nbm.2812
MRSCont = calculateGlobalConcentrations(MRSCont);
%% Atlas based analysis 
%   This function performs the atlas-based analysis of the MRSI data
[MRSCont] = OspreyMRSIAtlasAnalysis(MRSCont);

%% Save
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
MRSCont.flags.didOverview = 1;

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end


end