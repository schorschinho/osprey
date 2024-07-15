function [MRSCont] = OspreyOverview(MRSCont)
%% [MRSCont] = OspreyOverview(MRSCont)
%   This function creates te data structre needed for the overview panel
%   in the GUI. It sorts the data by groups and performs the needed
%   statistics.
%
%   USAGE:
%       MRSCont = OspreyOverview(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2019-02-19)
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
outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyOverview');


%%% 2. INITIALIZE VARIABLES %%%
%Getting the names of the SubSpectra and Fits


OrderNames = {'metab'};
OrderNamesFit = {'metab'};
if MRSCont.flags.hasMM
    OrderNames = horzcat(OrderNames, 'mm');
    OrderNamesFit = horzcat(OrderNamesFit, 'mm');
end

if MRSCont.flags.hasRef
    OrderNames = horzcat(OrderNames, 'ref');
    if ~strcmp(MRSCont.opts.fit.method,'LCModel')
        OrderNamesFit = horzcat(OrderNamesFit, 'ref');
    end
end

if MRSCont.flags.hasMMRef
    OrderNames = horzcat(OrderNames, 'mm_ref');
end

if MRSCont.flags.hasWater
    OrderNames = horzcat(OrderNames, 'w');
    if ~strcmp(MRSCont.opts.fit.method,'LCModel')
        OrderNamesFit = horzcat(OrderNamesFit, 'w');
    end
end
S = orderfields(MRSCont.processed,OrderNames);
dataPlotNames = fieldnames(S)';
NoSpec = length(fieldnames(MRSCont.processed));


for ss = 1:NoSpec
   SubSpecNamesStruct.(dataPlotNames{ss}) =  MRSCont.processed.(dataPlotNames{ss}){1}.names;
end
MRSCont.overview.SubSpecNamesStruct = SubSpecNamesStruct;

if MRSCont.flags.didFit
    S = orderfields(MRSCont.fit.results,OrderNamesFit);
    FitSpecNames = fieldnames(S)';
    NoFitSpecNames = length(FitSpecNames);
    if ~strcmp(MRSCont.opts.fit.method,'LCModel')
        for sf = 1 : NoFitSpecNames
            for sn = 1 : size(MRSCont.fit.results.(FitSpecNames{sf}).fitParams,3)
                for sb = 1 : size(MRSCont.fit.results.(FitSpecNames{sf}).fitParams,1)
                    if ~(strcmp(FitSpecNames{sf},'ref') ||strcmp(FitSpecNames{sf},'w'))
                        if ~isempty(MRSCont.fit.results.(FitSpecNames{sf}).fitParams{sb,1,sn})
                            FitSpecNamesStruct.(FitSpecNames{sf}){sb,sn} = MRSCont.fit.resBasisSet.(FitSpecNames{sf}).(MRSCont.info.(FitSpecNames{sf}).unique_ndatapoint_spectralwidth{1}){1,1,sn}.names{1};
                        end
                    else
                        FitSpecNamesStruct.(FitSpecNames{sf}){1} = MRSCont.fit.resBasisSet.(FitSpecNames{sf}).(MRSCont.info.(FitSpecNames{sf}).unique_ndatapoint_spectralwidth{1}){1,1,1}.names{1};
                    end
                end
            end
        end
    else
        FitSpecNamesStruct.(FitSpecNames{1}){1} = 'A';
        if MRSCont.flags.isMEGA
            FitSpecNamesStruct.(FitSpecNames{1}){2} = 'diff1';
        end
    end
    MRSCont.overview.FitSpecNamesStruct = FitSpecNamesStruct;
end



% shift = 0;


%%% 3. INTERPOLATION & NORMALIZATION %%%
% Starting with the processed data
OverviewTime = tic;

%Progress text for the GUI
if MRSCont.flags.isGUI && isfield(MRSCont.flags,'inProgress')
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end


%Interpolating spectra if needed to allow the calculation of mean and SD
%spectra
if (MRSCont.flags.isPRIAM == 1) && isfield(MRSCont.flags,'isPRIAM')
    Voxels = 2;
    for rr = 1 : Voxels
        for ss = 1 : NoSpec % Loop over Subspec
            for kk = 1 : MRSCont.nDatasets
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(dataPlotNames{ss}){kk},rr);
            end
        end
    end
else
    Voxels = 1;
    MRSCont.overview.Osprey.all_data_voxel_1 = MRSCont.processed;
end


fprintf('\n');
fprintf('Gathering spectra from subspectrum %d out of %d total subspectra...', 1, NoSpec);
for rr = 1 : Voxels
    for ss = 1 : NoSpec % Loop over Subspec
        msg = sprintf('Gathering spectra from subspectrum %d out of %d total subspectra...', ss, NoSpec);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
        if MRSCont.flags.isGUI && isfield(progressText,'String')
            set(progressText,'String' ,sprintf('Gathering spectra from subspectrum %d out of %d total subspectra...\n', ss, NoSpec));
            drawnow
        end
        for kk = 1 : MRSCont.nDatasets
            if MRSCont.processed.(dataPlotNames{ss}){1,kk}.sz(1) < MRSCont.info.(dataPlotNames{ss}).max_ndatapoint
                ppmRangeData        = MRSCont.processed.(dataPlotNames{ss}){1,MRSCont.info.(dataPlotNames{ss}).max_ndatapoint_ind}.ppm';
                ppmRangeDataToInt       = MRSCont.processed.(dataPlotNames{ss}){1,kk}.ppm;
                ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                if sum(ppmIsInDataRange) == 0
                    ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                end
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.specs      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.specs(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.ppm = ppmRangeData;
                if mod(size(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.dims.t),2)==0
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.fids=ifft(fftshift(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.dims.t),[],MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.dims.t);
                else
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.fids=ifft(circshift(fftshift(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.dims.t),1),[],MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){1,kk}.dims.t);
                end

            end
        end
    end
end

fprintf('\n... done.\n');
if MRSCont.flags.isGUI && isfield(progressText,'String')
    set(progressText,'String' ,sprintf('... done.'));
    pause(1);
end


if MRSCont.flags.didFit
    % Apply the same stpes to the fits
    fprintf('Gathering fit models from fit %d out of %d total fits...', 1, NoFitSpecNames);
    for rr = 1 : Voxels%Loop over
       for ss = 1 :NoFitSpecNames %Loop over fitted supsctra
       if strcmp(FitSpecNames{ss}, 'ref') || strcmp(FitSpecNames{ss}, 'w') % Water model
            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}) = cell(1,MRSCont.nDatasets(1),size(FitSpecNamesStruct.(FitSpecNames{ss}),2));
        else
            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}) = cell(size(FitSpecNamesStruct.(FitSpecNames{ss}),1),MRSCont.nDatasets(1),size(FitSpecNamesStruct.(FitSpecNames{ss}),2));
        end
        for sf = 1 : size(FitSpecNamesStruct.(FitSpecNames{ss}),2) %Loop over all fits
            for bf = 1 : size(FitSpecNamesStruct.(FitSpecNames{ss}),1) %Loop over all basis sets
                msg = sprintf('Gathering fit models from fit %d out of %d total fits...', ss, NoFitSpecNames);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                fprintf([reverseStr, msg]);
                if MRSCont.flags.isGUI && isfield(progressText,'String')
                    set(progressText,'String' ,sprintf('Gathering fit models from fit %d out of %d total fits...\n', ss, NoFitSpecNames));
                    drawnow
                end
                if ~isempty(FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf})
                    for kk = 1 : MRSCont.nDatasets(1) %Loop over all datasets
                        switch MRSCont.opts.fit.method %Which model was used
                        case 'Osprey'
                            if strcmp(FitSpecNames{ss}, 'ref') || strcmp(FitSpecNames{ss}, 'w') % Water model
                                % if water, use the water model
                                fitRangePPM = MRSCont.opts.fit.rangeWater;
                                if Voxels < 2
                                    dataToPlot  = MRSCont.processed.(FitSpecNames{ss}){kk};
                                    basisSet    = MRSCont.fit.resBasisSet.(FitSpecNames{ss}).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){1};
                                    % Get the fit parameters
                                    fitParams   = MRSCont.fit.results.(FitSpecNames{ss}).fitParams{kk};
                                else
                                    dataToPlot  = op_takeVoxel(MRSCont.processed.(FitSpecNames{ss}){kk},rr);
                                    basisSet    = MRSCont.fit.resBasisSet{rr}.(FitSpecNames{bf,sf}).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){1};
                                    % Get the fit parameters
                                    fitParams   = MRSCont.fit.results{rr}.(FitSpecNames{ss}).fitParams{kk};
                                end
                                % Pack up into structs to feed into the reconstruction functions
                                inputData.dataToFit                 = dataToPlot;
                                inputData.basisSet                  = basisSet;
                                if Voxels < 2
                                    inputSettings.scale                 = MRSCont.fit.scale{kk};
                                else
                                    inputSettings.scale                 = MRSCont.fit.scale{kk};
                                end
                                inputSettings.fitRangePPM           = fitRangePPM;
                                inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
                                % If water, extract and apply nonlinear parameters
                                [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fit      = ModelOutput.completeFit;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.ppm      = ModelOutput.ppm;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.data      = ModelOutput.data;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.res      = ModelOutput.residual;
                            else % if metabolite or MM data, use the metabolite model
                                fitRangePPM = MRSCont.opts.fit.range;
                                if strcmp(FitSpecNames{ss}, 'mm')
                                    fitRangePPM = [0.2 4.2];
                                end
                                if Voxels < 2

                                    dataToPlot  = op_takesubspec(MRSCont.processed.(FitSpecNames{ss}){kk},find(strcmp(MRSCont.processed.(FitSpecNames{ss}){kk}.names,FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf})));
                                    basisSet    = MRSCont.fit.resBasisSet.(FitSpecNames{ss}).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){bf,sf};
                                    if bf == 2 % We need to insert the subject specific MM basis function into the basis set
                                        if sf==1
                                            index = find(strcmp(MRSCont.processed.mm{kk}.names,'A_spline')); 
                                        end
                                        if sf==2
                                            index = find(strcmp(MRSCont.processed.mm{kk}.names,'diff1_spline')); 
                                        end
                                        mm_clean_spline = op_takesubspec(MRSCont.processed.mm{kk},index);
                                        mm_clean_spline               = op_zeropad(mm_clean_spline, 2);  
                                        ind = find(strcmp(basisSet.name,'MMExp'));
                                        basisSetfactor = op_freqrange(basisSet,0,1.2);
                                        mm_clean_spline_factor = op_freqrange(mm_clean_spline,0.7,1.1);
                                        factor = (max(real(basisSetfactor.specs(:,ind)))/max(real(mm_clean_spline_factor.specs)));
                                        mm_clean_spline = op_ampScale(mm_clean_spline,factor);
                                        basisSet.fids(:,ind) = mm_clean_spline.fids;
                                        basisSet.specs(:,ind) = mm_clean_spline.specs;
                                    end
                                    fitParams   = MRSCont.fit.results.(FitSpecNames{ss}).fitParams{bf,kk,sf};
                                else
                                   dataToPlot  = op_takeVoxel(MRSCont.processed.(dataPlotNames{ss}){kk},rr);
                                   fitParams   = MRSCont.fit.results{rr}.(FitSpecNames{ss}).fitParams{bf,kk,sf};
                                   basisSet    = MRSCont.fit.resBasisSet{rr}.(FitSpecNames{ss}).(['np_sw_' num2str(round(dataToPlot.sz(1))) '_' num2str(round(dataToPlot.spectralwidth))]){bf,sf};
                                end
                                % Pack up into structs to feed into the reconstruction functions
                                inputData.dataToFit                 = dataToPlot;
                                inputData.basisSet                  = basisSet;
                                if Voxels < 2
                                    inputSettings.scale                 = MRSCont.fit.scale{kk};
                                else
                                    inputSettings.scale                 = MRSCont.fit.scale{kk};
                                end
                                inputSettings.fitRangePPM           = fitRangePPM;
                                inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
                                inputSettings.fitStyle              = MRSCont.opts.fit.style;
                                inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
                                inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
                                inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
                                inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
                                inputSettings.concatenated.Subspec  = FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf};
                                if isfield(MRSCont.opts.fit,'GAP')
                                    inputSettings.GAP = MRSCont.opts.fit.GAP.(FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf});
                                else
                                    inputSettings.GAP = [];
                                end
                                if strcmp(inputSettings.fitStyle,'Concatenated')
                                    [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
                                else
                                    [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
                                end
                                if ~isnan(ModelOutput.completeFit(1)) %If the fit was succesfull
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fit      = ModelOutput.completeFit;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.baseline      = ModelOutput.baseline;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.ppm      =  ModelOutput.ppm;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.res      = ModelOutput.residual;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.data      = ModelOutput.data;
                                    if strcmp(FitSpecNames{ss}, 'mm') %re_mm loop over basis functions
                                        for n = 1 : (basisSet.nMets + basisSet.nMM)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.(['fit' basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                        end
                                        % tMM = all MM functions
                                        if basisSet.nMM > 0
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM  = sum(ModelOutput.indivMets(:,basisSet.nMets+1:end),2);
                                        else
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM =nan;
                                        end
                                        %section to write out MM_clean spectra
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.MM_clean = ModelOutput.data -sum(ModelOutput.indivMets(:,1:basisSet.nMets),2);
                                    else%re_mm
                                        for n = 1 : size(ModelOutput.indivMets,2) % loop over basis functions
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.(['fit' basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                        end
                                        % Add basis functions of metabolite combinations
                                        % tNAA = NAA + NAAG
                                        idx_NAA  = find(strcmp(basisSet.name,'NAA'));
                                        idx_NAAG  = find(strcmp(basisSet.name,'NAAG'));
                                        if isempty(idx_NAA) && isempty(idx_NAAG)
                                            % do nothing
                                        elseif isempty(idx_NAA) && ~isempty(idx_NAAG)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittNAA  = ModelOutput.indivMets(:,idx_NAAG);
                                        elseif ~isempty(idx_NAA) && isempty(idx_NAAG)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittNAA  = ModelOutput.indivMets(:,idx_NAA);
                                        elseif ~isempty(idx_NAA) && ~isempty(idx_NAAG)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittNAA  = ModelOutput.indivMets(:,idx_NAA) + ModelOutput.indivMets(:,idx_NAAG);
                                        end


                                        % tCr = Cr + tCr - CrCH2
                                        idx_Cr  = find(strcmp(basisSet.name,'Cr'));
                                        idx_PCr  = find(strcmp(basisSet.name,'PCr'));
                                        if isempty(idx_Cr) && isempty(idx_PCr)
                                            % do nothing
                                        elseif isempty(idx_Cr) && ~isempty(idx_PCr)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr  = ModelOutput.indivMets(:,idx_PCr);
                                        elseif ~isempty(idx_Cr) && isempty(idx_PCr)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr  = ModelOutput.indivMets(:,idx_Cr);
                                        elseif ~isempty(idx_Cr) && ~isempty(idx_PCr)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_PCr);
                                        end

                                        % if present, add CrCH2 model
                                        idx_CrCH2  = find(strcmp(basisSet.name,'CrCH2'));
                                        if ~isempty(idx_CrCH2)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr  = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr + ModelOutput.indivMets(:,idx_CrCH2);
                                        else
                                            % do nothing
                                        end

                                        % tCho = GPC + PCh
                                        idx_GPC  = find(strcmp(basisSet.name,'GPC'));
                                        idx_PCh  = find(strcmp(basisSet.name,'PCh'));
                                        if isempty(idx_GPC) && isempty(idx_PCh)
                                            % do nothing
                                        elseif isempty(idx_GPC) && ~isempty(idx_PCh)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCho  = ModelOutput.indivMets(:,idx_PCh);
                                        elseif ~isempty(idx_GPC) && isempty(idx_PCh)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCho  = ModelOutput.indivMets(:,idx_GPC);
                                        elseif ~isempty(idx_GPC) && ~isempty(idx_PCh)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCho  = ModelOutput.indivMets(:,idx_GPC) + ModelOutput.indivMets(:,idx_PCh);
                                        end

                                        % Glx = Glu + Gln
                                        idx_Glu  = find(strcmp(basisSet.name,'Glu'));
                                        idx_Gln  = find(strcmp(basisSet.name,'Gln'));
                                        if isempty(idx_Glu) && isempty(idx_Gln)
                                            % do nothing
                                        elseif isempty(idx_Glu) && ~isempty(idx_Gln)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fitGlx  = ModelOutput.indivMets(:,idx_Gln);
                                        elseif ~isempty(idx_Glu) && isempty(idx_Gln)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fitGlx  = ModelOutput.indivMets(:,idx_Glu);
                                        elseif ~isempty(idx_Glu) && ~isempty(idx_Gln)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fitGlx  = ModelOutput.indivMets(:,idx_Glu) + ModelOutput.indivMets(:,idx_Gln);
                                        end

                                        % tEA = PE + EA
                                        idx_PE  = find(strcmp(basisSet.name,'PE'));
                                        idx_EA  = find(strcmp(basisSet.name,'EA'));
                                        if isempty(idx_PE) && isempty(idx_Gln)
                                            % do nothing
                                        elseif isempty(idx_PE) && ~isempty(idx_EA)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittEA  = ModelOutput.indivMets(:,idx_EA);
                                        elseif ~isempty(idx_PE) && isempty(idx_EA)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittEA  = ModelOutput.indivMets(:,idx_PE);
                                        elseif ~isempty(idx_PE) && ~isempty(idx_EA)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittEA  = ModelOutput.indivMets(:,idx_PE) + ModelOutput.indivMets(:,idx_EA);
                                        end

                                        % tMM = all MM functions
                                        if MRSCont.opts.fit.fitMM == 1
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM  = sum(ModelOutput.indivMets(:,basisSet.nMets+1:end),2);
                                        else
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM =nan;
                                        end
                                    end %re_mm
                                else %if the fit was not succesful write nans into the corresponding fields
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fit      = nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.baseline      = nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.ppm      =  nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.res      = nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittNAA  = nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittCr  = nan;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.data      = nan;
                                end
                            end
                        case 'LCModel'
                                    if (MRSCont.flags.isPRIAM == 1)
                                        fitParams   = MRSCont.fit.results{rr}.(FitSpecNames{ss}).fitParams{bf,kk,sf};
                                    else
                                        fitParams   = MRSCont.fit.results.(FitSpecNames{ss}).fitParams{bf,kk,sf};
                                    end
                                    % Get the LCModel plots we previously extracted from .coord
                                    % etc.
                                    [ModelOutput] = fit_LCModelParamsToModel(fitParams);

                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fit      = ModelOutput.completeFit;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.baseline      = ModelOutput.baseline;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.ppm      =  ModelOutput.ppm';
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.res      = ModelOutput.residual;
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.data      = ModelOutput.data;
                                    for n = 1 : size(ModelOutput.indivMets,2) % loop over basis functions
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.(['fit' fitParams.name{n}])  = ModelOutput.indivMets(:,n);
                                    end
                                     % tMM = all MM functions
                                    if MRSCont.opts.fit.fitMM == 1
                                        %Find all MM or Lip functions that are not
                                        %combined
                                        idx_tMM = horzcat(find(contains(fitParams.name,'MM')), find(contains(fitParams.name,'Lip')));
                                        if ~isempty(idx_tMM)
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM  = sum(ModelOutput.indivMets(:,idx_tMM),2);
                                            else
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM =nan;
                                        end
                                    else
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(FitSpecNames{ss}){bf,kk,sf}.fittMM =nan;
                                    end
                        end
                    end
                end
            end
        end
       end
    end
    fprintf('\n... done.\n');
    if MRSCont.flags.isGUI  && isfield(progressText,'String')
        set(progressText,'String' ,sprintf('... done.'));
        pause(1);
    end

    ModelCombs = fieldnames(MRSCont.overview.Osprey.all_models_voxel_1);
    NoModelCombs = length(ModelCombs);
    for rr = 1 : Voxels
        for sc = 1 : NoModelCombs % Loop over all model combinations
            for sf = 1 : size(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}),3)
                for bf = 1 : size(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}),1)                                                                                                                                                                                                                                                           for bf = 1 : size(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}),1)
                        if isstruct(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf})
                            for kk = 1 : MRSCont.nDatasets
                                temp_fit_sz.(ModelCombs{sc})(bf,kk,sf)= length(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit);
                            end
                        end
                    end
                end
            end
            [max_point_fit.(ModelCombs{sc}),max_ind_fit.(ModelCombs{sc})] = max(temp_fit_sz.(ModelCombs{sc})(1,:,1));
        end
    end

    %Interpolating models if needed to allow the calculation of mean and SD
    %models
    fprintf('Interpolating fit models from fit %d out of %d total fits...', 1, NoModelCombs);
    for rr = 1 : Voxels
        for sc = 1 : NoModelCombs % Loop over all model combinations
           for sf = 1 : size(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}),3)
                for bf = 1 : size(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}),1)

                    msg = sprintf('Interpolating fit models from fit %d out of %d total fits...', sc, NoModelCombs);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    fprintf([reverseStr, msg]);
                    if MRSCont.flags.isGUI && isfield(progressText,'String')
                        set(progressText,'String' ,sprintf('Interpolating fit models from fit %d out of %d total fits...\n', sc, NoModelCombs));
                        drawnow
                    end

                    for kk = 1 : MRSCont.nDatasets %loop over all datasets
                        if isstruct(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf})
                            if length(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit) < max_point_fit.(ModelCombs{sc})
                                        ppmRangeData        = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){1,max_ind_fit.(ModelCombs{1})}.ppm';
                                        ppmRangeDataToInt       = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.ppm;
                                        ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                                        if sum(ppmIsInDataRange) == 0
                                            ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                                        end
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.data      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.data(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                        if ~(strcmp(ModelCombs{sc}, 'ref') || strcmp(ModelCombs{sc}, 'w'))
                                             MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.baseline = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.baseline(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                             names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf});
                                             for f = 6 : length(names)
                                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.(names{f})= interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.(names{f})(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                             end
                                        end
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.ppm = ppmRangeData';
                                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.res = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.data-MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit;
                                        if ~isempty(MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf}))
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.fit(ppmRangeData>MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(1) & ppmRangeData<MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(2)) = nan;
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.data(ppmRangeData>MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(1) & ppmRangeData<MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(2)) = nan;
                                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.baseline(ppmRangeData>MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(1) & ppmRangeData<MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(2)) = nan;
                                            names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf});
                                             for f = 6 : length(names)
                                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.(names{f})(ppmRangeData>MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(1) & ppmRangeData<MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(2)) = nan;
                                             end
                                             MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc}){bf,kk,sf}.res(ppmRangeData>MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(1) & ppmRangeData<MRSCont.opts.fit.GAP.(FitSpecNamesStruct.metab{sf})(2)) = nan;
                                        end
                            end
                        end
                    end
                end
            end
        end
    end
end

fprintf('\n... done.\n');
if MRSCont.flags.isGUI  && isfield(progressText,'String')
    set(progressText,'String' ,sprintf('... done.'));
    pause(1);
end

%%% 3. SCALING DATA  %%%
%Normalizing the data according to the scale value of the fit and normalize
%the models according to the tCr/tNAA amplitudes
fprintf('\nScaling data from dataset %d out of %d total datasets...', 1, MRSCont.nDatasets(1));
if ~MRSCont.flags.didFit
    [MRSCont] = osp_fitInitialise(MRSCont);
end

for rr = 1 : Voxels
    for ss = 1 : NoSpec % Loop over Subspec
        for kk = 1 : MRSCont.nDatasets
            if Voxels < 2
                scale                 = MRSCont.fit.scale{kk};
            else
                scale                 = MRSCont.fit.scale{kk};
            end
            msg = sprintf('Scaling data from dataset %d out of %d total datasets...', kk, MRSCont.nDatasets(1));
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            fprintf([reverseStr, msg]);
            if MRSCont.flags.isGUI && isfield(progressText,'String')
                set(progressText,'String' ,sprintf('Scaling data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets(1)));
                drawnow
            end
            if ~((strcmp(dataPlotNames{ss},'ref')||strcmp(dataPlotNames{ss},'mm_ref')||strcmp(dataPlotNames{ss},'w')) && strcmp(MRSCont.opts.fit.method, 'LCModel'))
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){kk}.specs = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){kk}.specs/scale;
            end
        end
    end
end
fprintf('\n... done.\n');
if MRSCont.flags.isGUI  && isfield(progressText,'String')
    set(progressText,'String' ,sprintf('... done.'));
    pause(1);
end


%%% 4. SORTING DATA  %%%
% Sort and group the data according to the stat.csv file. If no files is
% supplied a stat file with a single group is created. In addition, a grand
% mean is caclulated and the subject names are added into the stat csv file
% to allow an easier identification.

SepFileList = cell(1,length(MRSCont.files)); % Get all files
for kk = 1 : MRSCont.nDatasets(1)
    SepFileList{kk} =  split(MRSCont.files{kk}, filesep);
    ind = find(contains(lower(SepFileList{kk}),'sub'));
    if ~isempty(ind)
        subject{kk} = [SepFileList{kk}{ind(1)}]; % Create subject name list
    else
        subject{kk} = ['sub_' num2str(kk)];
    end

end

if MRSCont.flags.hasStatfile % Has stat file
    statFile = readtable(MRSCont.file_stat, 'Delimiter', ',','ReadVariableNames',1); % Load CSV input
    name = statFile.Properties.VariableNames;
    group_idx = find(strcmp(name,'group'));
    if isempty(group_idx) % No group supplied so create grand mean only
        MRSCont.overview.groups = ones(size(MRSCont.nDatasets,1),1);
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    else %Get grouping variable
        MRSCont.overview.groups = statFile{:,group_idx};
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    end
    if  ~strcmp(name,'subject') % No subject names stored in the container
        if length(subject)>1 && ~strcmp(subject{1},subject{2}) % Add names according to BIDS with subfolder
            if ~strcmp(name,'subject')
                statFile.subject = subject';
            end
        else
            if ~strcmp(name,'subject') % Add whole path as BIDS wasn't set up properly
                statFile.subject = MRSCont.files';
            end
        end

    end
else % No csv file supplied
    MRSCont.overview.groups = ones(MRSCont.nDatasets(1),1); %Create a single group
    MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    statFile = array2table(MRSCont.overview.groups,'VariableNames',{'group'});
    if length(subject)>1 && ~strcmp(subject{1},subject{2}) %Add names to the csv file
        statFile.subject = subject';
    else
        statFile.subject = MRSCont.files';
    end
    name = statFile.Properties.VariableNames;
end
if isfield(MRSCont, 'exclude') % If exclusions found in MRSCont, add to the table
    exclude = zeros(MRSCont.nDatasets(1),1);
    exclude(MRSCont.exclude) = 1;
    statFile.exclude = exclude;
end
statFile = addprop(statFile, {'VariableLongNames'}, {'variable'}); % add long name to table properties
% Loop over field names to populate descriptive fields of table for JSON export
for JJ = 1:length(name)
    switch name{JJ}
        case 'subject'
            statFile.Properties.CustomProperties.VariableLongNames{'subject'} = 'Subject'; %Write properties for json
            statFile.Properties.VariableDescriptions{'subject'} = 'Subject indetifier';
            statFile.Properties.VariableUnits{'subject'} = 'arbitrary';
        case 'group'
            statFile.Properties.CustomProperties.VariableLongNames{'group'} = 'Group'; %Write properties for json
            statFile.Properties.VariableDescriptions{'group'} = 'Sub-group the subject belongs to';
            statFile.Properties.VariableUnits{'group'} = 'arbitrary';
        case 'exclude'
            statFile.Properties.CustomProperties.VariableLongNames{'exclude'} = 'Excluded'; %Write properties for json
            statFile.Properties.VariableDescriptions{'exclude'} = 'Whether to exclude subject';
            statFile.Properties.VariableUnits{'exclude'} = 'arbitrary';
    end
end
osp_WriteBIDsTable(statFile,[MRSCont.outputFolder  filesep  'subject_names_and_excluded'])

%Exclude datasets based on the exclude field in the MRSConainer. THis can
%be triggered by pressing the left (remove) and right (add) arrow buttons
%in the listbox of the GUI
if isfield(MRSCont, 'exclude')
    if~isempty(MRSCont.exclude)
        MRSCont.overview.groups(MRSCont.exclude) = [];
        max_g = max(MRSCont.overview.groups);
        for kk = 1 : max(MRSCont.overview.groups)
            remove = 0;
            [~,idx] = find(MRSCont.overview.groups==kk);
            if isempty(idx) && kk < max_g
                for ll = 1 : length(MRSCont.overview.groups)
                    if MRSCont.overview.groups(ll) > kk
                        MRSCont.overview.groups(ll) = MRSCont.overview.groups(ll) -1;
                        if remove == 0
                            max_g = max_g -1;
                            remove =1;
                        end
                    end
                end
            end
        end
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    end
end

% Set up the different group names
MRSCont.overview.groupNames = cell(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    MRSCont.overview.groupNames{g} = ['Group ' num2str(g)];
end

% Sort the spectra according to the groups
for rr = 1 : Voxels
    for ss = 1 : NoSpec % Loop over subspectra
        for g = 1 : MRSCont.overview.NoGroups % loop over groups
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(['g_' num2str(g)]).(dataPlotNames{ss}) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss})(1,MRSCont.overview.groups == g);
        end
        MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).GMean.(dataPlotNames{ss}) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss})(1,MRSCont.overview.groups > 0);
    end
end

if MRSCont.flags.didFit
    for rr = 1 : Voxels
        % Sort the models according to the groups
        for sc = 1 : NoModelCombs % Loop over all model combinations
            for g = 1 : MRSCont.overview.NoGroups % loop over groups
                MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(['g_' num2str(g)]).(ModelCombs{sc}) = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc})(:,MRSCont.overview.groups == g,:);
            end
            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).GMean.(ModelCombs{sc}) = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).(ModelCombs{sc})(:,MRSCont.overview.groups > 0,:);
        end
    end
end

%%% 5. READ CORRELATION DATA INTO THE STRUCT %%%
% Open the stat file and add correlation measures to the MRSContainer
if MRSCont.flags.hasStatfile
    for cor = 1 : length(name)
        MRSCont.overview.corr.Names{cor} = name{cor};
    end
    for cor = 1 : length(name)
        MRSCont.overview.corr.Meas{cor} = statFile{:,cor};
        if isfield(MRSCont, 'exclude') % Exclude measures
            if~isempty(MRSCont.exclude)
                MRSCont.overview.corr.Meas{cor}(MRSCont.exclude) = [];
            end
        end
    end
end

%%% 6. CALCULATE MEAN AND SD SPECTRA FOR VISUALIZATION %%%
%Here we calculate the mean and SD spectra and fits for the overview plots

%Start with the spectra
for rr = 1 : Voxels
    for ss = 1 : NoSpec %loop over subspectra
        names = fields(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]));
        for g = 1 : length(names) % loop over groups
            tempSubSpec = zeros(length(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss})),MRSCont.info.(dataPlotNames{ss}).max_ndatapoint,MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss}){1, 1}.subspecs);
            if isempty(tempSubSpec)
                tempSubSpec = zeros(length(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss})),MRSCont.info.(dataPlotNames{ss}).max_ndatapoint,MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss}){1, 1}.rawSubspecs);
            end
            for kk = 1 : length(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss})) % Loop over datasets to generate a matrix
                try
                    tempSubSpec(kk,:,:) = MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(dataPlotNames{ss}){1,kk}.specs;
                catch
                    tempSubSpec(kk,:,:) = ones(1,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{1}){1,1}.sz(1)) *nan;
                end
            end
            %Calculate mean and SD
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(['mean_' dataPlotNames{ss}]) = squeeze(nanmean(real(tempSubSpec),1));
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(['sd_' dataPlotNames{ss}]) = squeeze(nanstd(real(tempSubSpec),1));
        end
        %Store ppm
        MRSCont.overview.Osprey.(['ppm_data_' dataPlotNames{ss}]) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(dataPlotNames{ss}){MRSCont.info.(dataPlotNames{ss}).max_ndatapoint_ind}.ppm;
    end
end

if MRSCont.flags.didFit
    %Do the same for the models
    for rr = 1 : Voxels
        for sc = 1 : NoModelCombs % Loop over all model combinations
            names = fields(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]));
            for g = 1 : length(names) %Loop over groups
               for sf = 1 : size(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}),3)
                for bf = 1 : size(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}),1)
                    if isstruct(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,1,sf})
                            tempSubSpec = zeros(size(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}),2),length(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,1,sf}.ppm));
                            tempSubRes = tempSubSpec;
                            tempSubdata = tempSubSpec;
                            tempSubBaseline = tempSubSpec;
                            tempInidivMetab = [];
                            for kk = 1 : size(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}),2) % Loop over datasets to generate a matrices
                              tempSubSpec(kk,:) = MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}.fit; %Fits
                              tempSubRes(kk,:) = MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}.res; % Residuals
                              tempSubdata(kk,:) = MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}.data; % spectra
                              if ~(strcmp(ModelCombs{sc}, 'ref') || strcmp(ModelCombs{sc}, 'w')) %Is not water
                                tempSubBaseline(kk,:) = MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}.baseline; % Baseline
                                fits = fields(MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}); % names of the basis functions
                                 for f = 6 : length(fits) % loop over basis functions
                                        tempInidivMetab.(fits{f})(kk,:)= MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,kk,sf}.(fits{f});
                                 end
                              end
                            end
                            %Calculate mean and SD
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['mean_fit_' ModelCombs{sc}])(bf,:,sf) = nanmean(real(tempSubSpec),1);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['sd_fit_' ModelCombs{sc}])(bf,:,sf) = nanstd(real(tempSubSpec),1);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['mean_res_' ModelCombs{sc}])(bf,:,sf) = nanmean(real(tempSubRes),1);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['sd_res_' ModelCombs{sc}])(bf,:,sf) = nanstd(real(tempSubRes),1);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['mean_data_' ModelCombs{sc}])(bf,:,sf) = nanmean(real(tempSubdata),1);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['sd_data_' ModelCombs{sc}])(bf,:,sf) = nanstd(real(tempSubdata),1);

                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['fit_' ModelCombs{sc}])(1:size(tempSubSpec,1),:,bf,sf) = real(tempSubSpec);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['res_' ModelCombs{sc}])(1:size(tempSubSpec,1),:,bf,sf) = real(tempSubRes);
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['data_' ModelCombs{sc}])(1:size(tempSubSpec,1),:,bf,sf) = real(tempSubdata);

                            if ~(strcmp(ModelCombs{sc}, 'ref') || strcmp(ModelCombs{sc}, 'w')) %Is not water
                                MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['mean_baseline_' ModelCombs{sc}])(bf,:,sf) = nanmean(real(tempSubBaseline),1);
                                MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['sd_baseline_' ModelCombs{sc}])(bf,:,sf) = nanstd(real(tempSubBaseline),1);
                                for f = 6 : length(fits) % loop over basis functions
                                        MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['mean_' fits{f} '_' ModelCombs{sc}])(bf,:,sf) = nanmean(real(tempInidivMetab.(fits{f})),1);
                                        MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['sd_' fits{f} '_' ModelCombs{sc}])(bf,:,sf) = nanstd(real(tempInidivMetab.(fits{f})),1);
                                end
                            end

                            %Store ppm
                            MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' ModelCombs{sc}])(bf,:,sf) = MRSCont.overview.Osprey.(['sort_models_voxel_' num2str(rr)]).(names{g}).(ModelCombs{sc}){bf,1,sf}.ppm;
                    end
                end
            end
           end
        end
    end
end

%%% 7. CLEAN UP AND SAVE %%%
% Set exit flags and version
MRSCont.flags.didOverview          = 1;
time = toc(OverviewTime);
fprintf('... done.\n Elapsed time %f seconds\n',time);
MRSCont.runtime.Overview = time;
fprintf('Runtime Breakdown................\n');
fprintf('OspreyLoad runtime: %f seconds\n',MRSCont.runtime.Load);
fprintf('OspreyProcess runtime: %f seconds\n',MRSCont.runtime.Proc);
if MRSCont.flags.didFit
    fprintf('OspreyFit runtime: %f seconds\n',MRSCont.runtime.Fit);
    fprintf('\tOspreyFit metab runtime: %f seconds\n',MRSCont.runtime.FitMet);
    if isfield(MRSCont.runtime, 'FitRef')
        fprintf('\tOspreyFit reference runtime: %f seconds\n',MRSCont.runtime.FitRef);
    end
    if isfield(MRSCont.runtime, 'FitWater')
        fprintf('\tOspreyFit water runtime: %f seconds\n',MRSCont.runtime.FitWater);
    end
else
    MRSCont.runtime.Fit = 0;
end
if ~MRSCont.flags.didQuantify
    MRSCont.runtime.Quantify = 0;
end
MRSCont.runtime.All = MRSCont.runtime.Load +MRSCont.runtime.Proc+MRSCont.runtime.Fit+MRSCont.runtime.Quantify+MRSCont.runtime.Overview;
if isfield(MRSCont.runtime, 'Coreg')
    MRSCont.runtime.All = MRSCont.runtime.All + MRSCont.runtime.Coreg;
    fprintf('OspreyCoreg runtime: %f seconds\n',MRSCont.runtime.Coreg);
end
if isfield(MRSCont.runtime, 'Seg')
    MRSCont.runtime.All = MRSCont.runtime.All + MRSCont.runtime.Seg;
    fprintf('OspreySeg runtime: %f seconds\n',MRSCont.runtime.Seg);
end
fprintf('OspreyOverview runtime: %f seconds\n',MRSCont.runtime.Overview);
fprintf('Full Osprey runtime: %f seconds\n',MRSCont.runtime.All);

diary off
% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;


% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyOverview')
end

% Create the MRSinMRS markdown
if MRSCont.flags.didFit && MRSCont.flags.didQuantify
    [MRSCont] = OspreyMinReport(MRSCont);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

if MRSCont.flags.isGUI  && isfield(progressText,'String')
    set(progressText,'String' ,sprintf('\n Elapsed time %f seconds',time));
    pause(1);
end

if MRSCont.opts.exportParams.flag==1
    osp_exportParams(MRSCont)
end

end
