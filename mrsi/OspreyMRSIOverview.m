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

MRSCont = calculateGlobalConcentrations(MRSCont);
%% Atlas based analysis 

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