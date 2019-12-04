function [MRSCont] = OspreyOverview(MRSCont)
%% [MRSCont] = OspreyOverview(MRSCont)
%   This function creates te data structre needed for the overview panel
%   in the GUI. It sorts the data by groups and performs the needed
%   statistics.
%
%   USAGE:
%       MRSCont = OspreyLoad(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2019-02-19)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-11-11: First version of the code.

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
    if nargin<1
        error('ERROR: no input Osprey container specified.  Aborting!!');
    end

%%% 2. INITIALIZE VARIABLES %%%
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));

%%% 3. ZEROFILLING & NORMALIZATION %%%
MRSCont.overview.all_data = MRSCont.processed;
temp_sz = zeros(1,MRSCont.nDatasets);
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
        temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
    end
    max_point = max(temp_sz);
    for kk = 1 : MRSCont.nDatasets
        if MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1) < max_point
            zf = max_point/MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}=op_zeropad(MRSCont.processed.(SubSpecNames{ss}){1,kk},zf);
        end
    end
end
MRSCont.overview.all_data_NAAnormalized = MRSCont.overview.all_data;
for kk = 1 : MRSCont.nDatasets
    if isfield(MRSCont, 'fit')
        if MRSCont.flags.isUnEdited
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/max(real(MRSCont.overview.all_data.A{1,kk}.specs(MRSCont.overview.all_data.A{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.A{1,kk}.ppm > 2.1)));
                end
            end
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.A{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.A{1,kk}.ppm > 2.1)));
        end
        if MRSCont.flags.isMEGA
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/max(real(MRSCont.overview.all_data.B{1,kk}.specs(MRSCont.overview.all_data.B{1,kk}.ppm > 1.9 && MRSCont.overview.all_data.B{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/max(real(MRSCont.overview.all_data.B{1,kk}.specs(MRSCont.overview.all_data.B{1,kk}.ppm > 1.9 && MRSCont.overview.all_data.B{1,kk}.ppm > 2.1)));
                end
            end
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.B{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.B{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.B{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.B{1,kk}.ppm > 2.1)));
        end
        if MRSCont.flags.isHERMES
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                end
            end
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
         end
        if MRSCont.flags.isHERCULES
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/max(real(MRSCont.overview.all_data.C{1,kk}.specs(MRSCont.overview.all_data.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data.C{1,kk}.ppm > 2.1)));
                end
            end
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
            MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs/max(real(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs(MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 1.9 & MRSCont.overview.all_data_NAAnormalized.C{1,kk}.ppm > 2.1)));
        end
    else
        error('This script works only on fully processed data. Run the whole Osprey pipeline first. Seg/Coreg is not needed')
    end
end

%%% 4. SORTING DATA  %%%
if MRSCont.flags.hasStatfile
    statCSV = readtable(MRSCont.file_stat, 'Delimiter', ',');
    group_idx = find(strcmp(statCSV{:,end},'group'));
    if isempty(group_idx)
        MRSCont.overview.groups = ones(MRSCont.nDatasets,1);
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    else
        MRSCont.overview.groups = statCSV{:,group_idx};
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    end
else
    MRSCont.overview.groups = ones(MRSCont.nDatasets,1);
    MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
end
MRSCont.overview.groupNames = cell(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    MRSCont.overview.groupNames{g} = ['Group ' num2str(g)];
end

for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}) = MRSCont.overview.all_data.(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}) = MRSCont.overview.all_data_NAAnormalized.(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
    end
end

%%% 5. CALCULATE MEAN AND SD SPECTRA FOR VISUALIZATION %%%

for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        tempSubSpec = zeros(length(MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss})),MRSCont.overview.all_data.(SubSpecNames{1}){1,1}.sz(1));
        for kk = 1 : length(MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}))
          tempSubSpec(kk,:) = MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}){1,kk}.specs;
        end
        MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' SubSpecNames{ss}]) = mean(real(tempSubSpec),1);
        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' SubSpecNames{ss}]) = std(real(tempSubSpec),1);
    end
MRSCont.overview.(['ppm_' SubSpecNames{ss}]) = MRSCont.overview.all_data.(SubSpecNames{ss}){1,1}.ppm;
end

for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        tempSubSpec = zeros(length(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss})),MRSCont.overview.all_data_NAAnormalized.(SubSpecNames{1}){1,1}.sz(1));
        for kk = 1 : length(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}))
          tempSubSpec(kk,:) = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}){1,kk}.specs;
        end
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' SubSpecNames{ss}]) = mean(real(tempSubSpec),1);
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' SubSpecNames{ss}]) = std(real(tempSubSpec),1);
    end
end

%%% 6. READ CORRELATION DATA IN THE STRUCT %%%
if MRSCont.flags.hasStatfile
    cor = 1;
    while ~strcmp(statCSV{cor,end},'')
        MRSCont.overview.corr.Names{cor} = statCSV{cor,end};
        cor = cor + 1;
    end
    for cor = 1 : size(statCSV,2)-1
        MRSCont.overview.corr.Meas{cor} = statCSV{:,cor};
    end
end

%%% 7. CLEAN UP AND SAVE %%%
% Set exit flags
MRSCont.flags.didOverview          = 1;

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
save(fullfile(outputFolder, outputFile), 'MRSCont');

end
