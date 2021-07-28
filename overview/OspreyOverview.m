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
osp_CheckRunPreviousModule(MRSCont, 'OspreyOverview');
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check ('OspreyOverview',MRSCont.flags.isGUI);


%%% 2. INITIALIZE VARIABLES %%%
%Getting the names of the SubSpectra and Fits
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));

if MRSCont.flags.didFit
    if iscell(MRSCont.fit.results) % Is DualVoxel
        FitNames = fieldnames(MRSCont.fit.results{1});
        NoFit = length(fieldnames(MRSCont.fit.results{1}));
    else %Is SVS
        FitNames = fieldnames(MRSCont.fit.results);
        NoFit = length(fieldnames(MRSCont.fit.results));
    end

    dataPlotNames = FitNames;
    tempFitNames = FitNames;
end

shift = 0;

if MRSCont.flags.didFit
    %Getting the final model names (needed for concatenated fits)
    for sf = 1 : NoFit
        switch MRSCont.opts.fit.method
            case 'Osprey'
                switch FitNames{sf}
                    case 'off'
                        dataPlotNames{sf} = 'A';
                    case 'conc'
                        if MRSCont.flags.isMEGA
                            dataPlotNames{sf} = 'diff1';
                            dataPlotNames{sf+1} = 'sum';
                            tempFitNames{sf} = 'conc';
                            tempFitNames{sf+1} = 'conc';
                            shift = 1;
                        end
                        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                            dataPlotNames{sf} = 'diff1';
                            dataPlotNames{sf+1} = 'diff2';
                            dataPlotNames{sf+2} = 'sum';
                            tempFitNames{sf} = 'conc';
                            tempFitNames{sf+1} = 'conc';
                            tempFitNames{sf+2} = 'conc';
                            shift = 2;
                        end
                    otherwise
                        dataPlotNames{sf + shift} = FitNames{sf};
                        tempFitNames{sf + shift} = FitNames{sf};
                end
            case 'OspreyAsym'
                switch FitNames{sf}
                    case 'off'
                        dataPlotNames{sf} = 'A';
                    case 'conc'
                        if MRSCont.flags.isMEGA
                            dataPlotNames{sf} = 'diff1';
                            dataPlotNames{sf+1} = 'sum';
                            tempFitNames{sf} = 'conc';
                            tempFitNames{sf+1} = 'conc';
                            shift = 1;
                        end
                        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                            dataPlotNames{sf} = 'diff1';
                            dataPlotNames{sf+1} = 'diff2';
                            dataPlotNames{sf+2} = 'sum';
                            tempFitNames{sf} = 'conc';
                            tempFitNames{sf+1} = 'conc';
                            tempFitNames{sf+2} = 'conc';
                            shift = 2;
                        end
                    otherwise
                        dataPlotNames{sf + shift} = FitNames{sf};
                        tempFitNames{sf + shift} = FitNames{sf};
                end            
            case 'LCModel'
                switch FitNames{sf}
                    case 'off'
                        dataPlotNames{sf} = 'A';
                    otherwise
                        dataPlotNames{sf} = FitNames{sf};
                end
        end
    end
    FitNames = tempFitNames;
    NoFit = length(FitNames);
end

%%% 3. INTERPOLATION & NORMALIZATION %%%
% Starting with the processed data 
OverviewTime = tic;

%Progress text for the GUI
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end


%Interpolating spectra if needed to allow the calculation of mean and SD
%spectra
if (MRSCont.flags.isPRIAM == 1) && isfield(MRSCont.flags,'isPRIAM')
    Voxels = 2;
    for rr = 1 : Voxels
        for ss = 1 : NoSubSpec % Loop over Subspec
            for kk = 1 : MRSCont.nDatasets
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){kk} = op_takeVoxel(MRSCont.processed.(SubSpecNames{ss}){kk},rr);
            end
        end
    end
else
    Voxels = 1;
    MRSCont.overview.Osprey.all_data_voxel_1 = MRSCont.processed;
end
fprintf('\n');
fprintf('Gathering spectra from subspectrum %d out of %d total subspectra...', 1, NoSubSpec);
for rr = 1 : Voxels
    for ss = 1 : NoSubSpec % Loop over Subspec
        msg = sprintf('Gathering spectra from subspectrum %d out of %d total subspectra...', ss, NoSubSpec);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, msg]);
        if MRSCont.flags.isGUI && isfield(progressText,'String')
            set(progressText,'String' ,sprintf('Gathering spectra from subspectrum %d out of %d total subspectra...\n', ss, NoSubSpec));
            drawnow
        end
        for kk = 1 : MRSCont.nDatasets
            if MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1) < MRSCont.info.(SubSpecNames{ss}).max_ndatapoint
                ppmRangeData        = MRSCont.processed.(SubSpecNames{ss}){1,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind}.ppm';
                ppmRangeDataToInt       = MRSCont.processed.(SubSpecNames{ss}){1,kk}.ppm;
                ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                if sum(ppmIsInDataRange) == 0
                    ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                end
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm = ppmRangeData;
                if mod(size(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.dims.t),2)==0
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.fids=ifft(fftshift(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.dims.t),[],MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.dims.t);
                else
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.fids=ifft(circshift(fftshift(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.dims.t),1),[],MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.dims.t);
                end

            end
        end
    end
end


% Align the spectra according to the NAA peak
for rr = 1 : Voxels
    for ss = 1 : NoSubSpec
        for kk = 1 : MRSCont.nDatasets
         %Find the ppm of the maximum peak magnitude within the given range:
            if MRSCont.flags.isUnEdited
                %Find the ppm of the maximum peak magnitude within the given range:
                ppmindex=find(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm<2.1)));
                ppmrange=MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm<2.1);
                ppmmax=ppmrange(ppmindex);
                refShift=(ppmmax-2.013);
            end
            if MRSCont.flags.isMEGA
                ppmindex=find(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1)));
                ppmrange=MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1);
                ppmmax=ppmrange(ppmindex);
                refShift=(ppmmax-2.013);
            end
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                ppmindex=find(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1)));
                ppmrange=MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm(MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.ppm<2.1);
                ppmmax=ppmrange(ppmindex);
                refShift=(ppmmax-2.013);
            end

         MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){1,kk}.ppm - refShift;
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
    fprintf('Gathering fit models from fit %d out of %d total fits...', 1, NoFit);
    for rr = 1 : Voxels%Loop over
        for sf = 1 : NoFit %Loop over all fits
            msg = sprintf('Gathering fit models from fit %d out of %d total fits...', sf, NoFit);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            fprintf([reverseStr, msg]);
            if MRSCont.flags.isGUI && isfield(progressText,'String')
                set(progressText,'String' ,sprintf('Gathering fit models from fit %d out of %d total fits...\n', sf, NoFit));
                drawnow
            end
            for kk = 1 : MRSCont.nDatasets %Loop over all datasets
                switch MRSCont.opts.fit.method %Which model was used
                case 'Osprey'
                    if strcmp((FitNames{sf}), 'ref') || strcmp((FitNames{sf}), 'w') % Water model 
                        % if water, use the water model
                        fitRangePPM = MRSCont.opts.fit.rangeWater;
                        if Voxels < 2
                            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
                            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
                        else
                            dataToPlot  = op_takeVoxel(MRSCont.processed.(dataPlotNames{sf}){kk},rr);
                            basisSet    = MRSCont.fit.resBasisSet{rr}.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);                            
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results{rr}.(FitNames{sf}).fitParams{kk};
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
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      = ModelOutput.ppm;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;
                    else % if metabolite or MM data, use the metabolite model           
                        fitRangePPM = MRSCont.opts.fit.range;
                        if Voxels < 2
                            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
                            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);
                            fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};                            
                        else
                           dataToPlot  = op_takeVoxel(MRSCont.processed.(dataPlotNames{sf}){kk},rr);
                           fitParams   = MRSCont.fit.results{rr}.(FitNames{sf}).fitParams{kk}; 
                           basisSet    = MRSCont.fit.resBasisSet{rr}.(FitNames{sf}).(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);
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
                        inputSettings.concatenated.Subspec  = dataPlotNames{sf};
                        if strcmp(inputSettings.fitStyle,'Concatenated')
                            [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
                        else
                            [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
                        end
                        if ~isnan(ModelOutput.completeFit) %If the fit was succesful
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = ModelOutput.baseline;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  ModelOutput.ppm;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
                            if strcmp(FitNames{sf}, 'mm') %re_mm loop over basis functions
                                for n = 1 : 4 + MRSCont.fit.basisSet.nMM
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' MRSCont.fit.basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                end
                                idx_NAA  = 4;
                                idx_Cr  = 1;
                                idx_CrCH2  = 2;
                                if ~isempty(idx_CrCH2)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_CrCH2);
                                else
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) ;
                                end
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,5:end),2);
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitNAA  = sum(ModelOutput.indivMets(:,4),2);

                                end
                            %section to write out MM_clean spectra
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.MM_clean = ModelOutput.data -sum(ModelOutput.indivMets(:,1:4),2);
                            else%re_mm
                                for n = 1 : size(ModelOutput.indivMets,2) % loop over basis functions
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                end
                                % Add basis functions of metabolite combinations 
                                % tNAA = NAA + NAAG
                                idx_NAA  = find(strcmp(basisSet.name,'NAA'));
                                idx_NAAG  = find(strcmp(basisSet.name,'NAAG'));
                                if isempty(idx_NAA) && isempty(idx_NAAG)
                                    % do nothing
                                elseif isempty(idx_NAA) && ~isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAAG);
                                elseif ~isempty(idx_NAA) && isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAA);
                                elseif ~isempty(idx_NAA) && ~isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAA) + ModelOutput.indivMets(:,idx_NAAG);
                                end


                                % tCr = Cr + tCr - CrCH2
                                idx_Cr  = find(strcmp(basisSet.name,'Cr'));
                                idx_PCr  = find(strcmp(basisSet.name,'PCr'));
                                if isempty(idx_Cr) && isempty(idx_PCr)
                                    % do nothing
                                elseif isempty(idx_Cr) && ~isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_PCr);
                                elseif ~isempty(idx_Cr) && isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr);
                                elseif ~isempty(idx_Cr) && ~isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_PCr);
                                end

                                % if present, add CrCH2 model
                                idx_CrCH2  = find(strcmp(basisSet.name,'CrCH2'));
                                if ~isempty(idx_CrCH2)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr + ModelOutput.indivMets(:,idx_CrCH2);
                                else
                                    % do nothing
                                end

                                % tCho = GPC + PCh
                                idx_GPC  = find(strcmp(basisSet.name,'GPC'));
                                idx_PCh  = find(strcmp(basisSet.name,'PCh'));
                                if isempty(idx_GPC) && isempty(idx_PCh)
                                    % do nothing
                                elseif isempty(idx_GPC) && ~isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_PCh);
                                elseif ~isempty(idx_GPC) && isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_GPC);
                                elseif ~isempty(idx_GPC) && ~isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_GPC) + ModelOutput.indivMets(:,idx_PCh);
                                end

                                % Glx = Glu + Gln
                                idx_Glu  = find(strcmp(basisSet.name,'Glu'));
                                idx_Gln  = find(strcmp(basisSet.name,'Gln'));
                                if isempty(idx_Glu) && isempty(idx_Gln)
                                    % do nothing
                                elseif isempty(idx_Glu) && ~isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Gln);
                                elseif ~isempty(idx_Glu) && isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Glu);
                                elseif ~isempty(idx_Glu) && ~isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Glu) + ModelOutput.indivMets(:,idx_Gln);
                                end

                                % tEA = PE + EA
                                idx_PE  = find(strcmp(basisSet.name,'PE'));
                                idx_EA  = find(strcmp(basisSet.name,'EA'));
                                if isempty(idx_PE) && isempty(idx_Gln)
                                    % do nothing
                                elseif isempty(idx_PE) && ~isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_EA);
                                elseif ~isempty(idx_PE) && isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_PE);
                                elseif ~isempty(idx_PE) && ~isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_PE) + ModelOutput.indivMets(:,idx_EA);
                                end

                                % tMM = all MM functions
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,basisSet.nMets+1:end),2);
                                end
                            end %re_mm
                        else %if the fit was not succesful write nans into the corresponding fields
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = nan;
                        end
                    end
                case 'OspreyAsym'
                        if strcmp((FitNames{sf}), 'ref') || strcmp((FitNames{sf}), 'w') % Water model 
                        % if water, use the water model
                       if Voxels < 2
                            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
                            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
                        else
                            dataToPlot  = op_takeVoxel(MRSCont.processed.(dataPlotNames{sf}){kk},rr);
                            basisSet    = MRSCont.fit.resBasisSet{rr}.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);                            
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results{rr}.(FitNames{sf}).fitParams{kk};
                        end
                        % Get the fit parameters
                        fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
                        % Pack up into structs to feed into the reconstruction functions
                        inputData.dataToFit                 = dataToPlot;
                        inputData.basisSet                  = basisSet;
                        inputSettings.scale                 = MRSCont.fit.scale{kk};
                        inputSettings.fitRangePPM           = fitRangePPM;
                        inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
                        % If water, extract and apply nonlinear parameters
                        [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      = ModelOutput.ppm;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;
                    else % if metabolite or MM data, use the metabolite model           
                        fitRangePPM = MRSCont.opts.fit.range;
                        if Voxels < 2
                            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
                            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
                        else
                            dataToPlot  = op_takeVoxel(MRSCont.processed.(dataPlotNames{sf}){kk},rr);
                            basisSet    = MRSCont.fit.resBasisSet{rr}.(FitNames{sf}).water.(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]);                            
                            % Get the fit parameters
                            fitParams   = MRSCont.fit.results{rr}.(FitNames{sf}).fitParams{kk};
                        end
                        % Pack up into structs to feed into the reconstruction functions
                        inputData.dataToFit                 = dataToPlot;
                        inputData.basisSet                  = basisSet;
                        if Voxels < 2
                            inputSettings.scale                 = MRSCont.fit.scale{kk};
                        else
                            inputSettings.scale                 = MRSCont.fit.scale{kk}(rr);                        
                        end
                        inputSettings.fitRangePPM           = fitRangePPM;
                        inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
                        inputSettings.fitStyle              = MRSCont.opts.fit.style;
                        inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
                        inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
                        inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
                        inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
                        inputSettings.concatenated.Subspec  = dataPlotNames{sf};
                        if strcmp(inputSettings.fitStyle,'Concatenated')
                            [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
                        else
                            [ModelOutput] = fit_OspreyAsymParamsToModel(inputData, inputSettings, fitParams);
                        end
                        if ~isnan(ModelOutput.completeFit) %If the fit was succesful
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = ModelOutput.baseline;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  ModelOutput.ppm;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
                            if strcmp(FitNames{sf}, 'mm') %re_mm loop over basis functions
                                for n = 1 : 4 + MRSCont.fit.basisSet.nMM
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' MRSCont.fit.basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                end
                                idx_NAA  = 4;
                                idx_Cr  = 1;
                                idx_CrCH2  = 2;
                                if ~isempty(idx_CrCH2)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_CrCH2);
                                else
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) ;
                                end
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,5:end),2);
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitNAA  = sum(ModelOutput.indivMets(:,4),2);

                                end
                            %section to write out MM_clean spectra
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.MM_clean = ModelOutput.data -sum(ModelOutput.indivMets(:,1:4),2);
                            else%re_mm
                                for n = 1 : size(ModelOutput.indivMets,2) % loop over basis functions
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' basisSet.name{n}])  = ModelOutput.indivMets(:,n);
                                end
                                % Add basis functions of metabolite combinations 
                                % tNAA = NAA + NAAG
                                idx_NAA  = find(strcmp(basisSet.name,'NAA'));
                                idx_NAAG  = find(strcmp(basisSet.name,'NAAG'));
                                if isempty(idx_NAA) && isempty(idx_NAAG)
                                    % do nothing
                                elseif isempty(idx_NAA) && ~isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAAG);
                                elseif ~isempty(idx_NAA) && isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAA);
                                elseif ~isempty(idx_NAA) && ~isempty(idx_NAAG)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAA) + ModelOutput.indivMets(:,idx_NAAG);
                                end


                                % tCr = Cr + tCr - CrCH2
                                idx_Cr  = find(strcmp(basisSet.name,'Cr'));
                                idx_PCr  = find(strcmp(basisSet.name,'PCr'));
                                if isempty(idx_Cr) && isempty(idx_PCr)
                                    % do nothing
                                elseif isempty(idx_Cr) && ~isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_PCr);
                                elseif ~isempty(idx_Cr) && isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr);
                                elseif ~isempty(idx_Cr) && ~isempty(idx_PCr)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_PCr);
                                end

                                % if present, add CrCH2 model
                                idx_CrCH2  = find(strcmp(basisSet.name,'CrCH2'));
                                if ~isempty(idx_CrCH2)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr + ModelOutput.indivMets(:,idx_CrCH2);
                                else
                                    % do nothing
                                end

                                % tCho = GPC + PCh
                                idx_GPC  = find(strcmp(basisSet.name,'GPC'));
                                idx_PCh  = find(strcmp(basisSet.name,'PCh'));
                                if isempty(idx_GPC) && isempty(idx_PCh)
                                    % do nothing
                                elseif isempty(idx_GPC) && ~isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_PCh);
                                elseif ~isempty(idx_GPC) && isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_GPC);
                                elseif ~isempty(idx_GPC) && ~isempty(idx_PCh)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_GPC) + ModelOutput.indivMets(:,idx_PCh);
                                end

                                % Glx = Glu + Gln
                                idx_Glu  = find(strcmp(basisSet.name,'Glu'));
                                idx_Gln  = find(strcmp(basisSet.name,'Gln'));
                                if isempty(idx_Glu) && isempty(idx_Gln)
                                    % do nothing
                                elseif isempty(idx_Glu) && ~isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Gln);
                                elseif ~isempty(idx_Glu) && isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Glu);
                                elseif ~isempty(idx_Glu) && ~isempty(idx_Gln)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_Glu) + ModelOutput.indivMets(:,idx_Gln);
                                end

                                % tEA = PE + EA
                                idx_PE  = find(strcmp(basisSet.name,'PE'));
                                idx_EA  = find(strcmp(basisSet.name,'EA'));
                                if isempty(idx_PE) && isempty(idx_Gln)
                                    % do nothing
                                elseif isempty(idx_PE) && ~isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_EA);
                                elseif ~isempty(idx_PE) && isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_PE);
                                elseif ~isempty(idx_PE) && ~isempty(idx_EA)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittEA  = ModelOutput.indivMets(:,idx_PE) + ModelOutput.indivMets(:,idx_EA);
                                end

                                % tMM = all MM functions
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,basisSet.nMets+1:end),2);
                                end
                            end %re_mm
                        else %if the fit was not succesful write nans into the corresponding fields
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = nan;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = nan;
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

    for rr = 1 : Voxels
        for sf = 1 : NoFit
            for kk = 1 : MRSCont.nDatasets
                temp_fit_sz.([FitNames{sf} '_' dataPlotNames{sf}])(1,kk)= length(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit);
            end
            [max_point_fit.([FitNames{sf} '_' dataPlotNames{sf}]),max_ind_fit.([FitNames{sf} '_' dataPlotNames{sf}])] = max(temp_fit_sz.([FitNames{sf} '_' dataPlotNames{sf}]));
        end
    end

    %Interpolating models if needed to allow the calculation of mean and SD
    %models
    fprintf('Interpolating fit models from fit %d out of %d total fits...', 1, NoFit);
    for rr = 1 : Voxels
        for sf = 1 : NoFit % loop over all fits
            msg = sprintf('Interpolating fit models from fit %d out of %d total fits...', sf, NoFit);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            fprintf([reverseStr, msg]);
            if MRSCont.flags.isGUI && isfield(progressText,'String')
                set(progressText,'String' ,sprintf('Interpolating fit models from fit %d out of %d total fits...\n', sf, NoFit));
                drawnow
            end
            for kk = 1 : MRSCont.nDatasets %loop over all datasets
                if length(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit) < max_point_fit.([FitNames{sf} '_' dataPlotNames{sf}])
                            ppmRangeData        = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,max_ind_fit.([FitNames{sf} '_' dataPlotNames{sf}])}.ppm';
                            ppmRangeDataToInt       = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm;
                            ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                            if sum(ppmIsInDataRange) == 0
                                ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                            end
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                            if ~(strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w'))
                                 MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                 names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk});
                                 for f = 6 : length(names)
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(names{f})= interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(names{f})(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                 end
                            end
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm = ppmRangeData';
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data-MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit;
                end
            end
        end
    end

    % Align the spectra according to the NAA peak
    for rr = 1 : Voxels
        for sf = 1 : NoFit % loop over all fits
            if ~strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w')
                for kk = 1 : MRSCont.nDatasets % loop over all data sets
                    naa = 1;
                     %Find the ppm of the maximum peak magnitude within the given range:
                     if MRSCont.flags.isUnEdited
                        ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1)));
                        ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1);
                     end
                     if MRSCont.flags.isMEGA
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1);
                        else
                            ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<2.1);
                        end
                        if isempty(ppmindex)
                            ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1)));
                            ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1);
                            naa = 0;
                        end
                     end
                     if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.ppm<2.1);
                        else
                            ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.ppm<2.1);
                        end
                     end
                    ppmmax=ppmrange(ppmindex);
                    if naa == 1
                        MRSCont.overview.Osprey.refShift(kk)=(ppmmax-2.013); %ref shift value
                    else
                        MRSCont.overview.Osprey.refShift(kk)=(ppmmax-3.03); %ref shift value
                    end
                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm - MRSCont.overview.Osprey.refShift(kk);
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
fprintf('\nScaling data from dataset %d out of %d total datasets...', 1, MRSCont.nDatasets);
if MRSCont.flags.didFit
    for rr = 1 : Voxels
        for kk = 1 : MRSCont.nDatasets
            if Voxels < 2
                scale                 = MRSCont.fit.scale{kk};
            else
                scale                 = MRSCont.fit.scale{kk};                        
            end
            msg = sprintf('Scaling data from dataset %d out of %d total datasets...', kk, MRSCont.nDatasets);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            fprintf([reverseStr, msg]);
            if MRSCont.flags.isGUI && isfield(progressText,'String')
                set(progressText,'String' ,sprintf('Scaling data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
                drawnow
            end
            if isfield(MRSCont, 'quantify')
                if MRSCont.flags.isUnEdited 
                    if MRSCont.flags.hasMM %re_mm
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).mm{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).mm{1,kk}.specs/scale; %re_mm
                        names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).mm_mm{1,kk});
                        for f = 1 : length(names)
                             if ~strcmp(names{f},'ppm')
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).mm_mm{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).mm_mm{1,kk}.(names{f})/(MRSCont.fit.results.off.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.off.fitParams{1,kk}.ampl(idx_Cr));
                             end
                        end
                    end %re_mm
                    if MRSCont.flags.hasRef
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit;
                    end
                    if MRSCont.flags.hasWater
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit;
                    end
                    names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk});
                     for f = 1 : length(names)
                         if ~strcmp(names{f},'ppm')
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.(names{f});
                         end
                     end
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                end

                if MRSCont.flags.isMEGA
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk});
                             for f = 1 : length(names)
                                 if ~strcmp(names{f},'ppm')
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.(names{f});
                                 end
                             end
                              names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk});
                             for f = 1 : length(names)
                                 if ~strcmp(names{f},'ppm')
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.(names{f});
                                 end
                             end
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs/scale;
                            if MRSCont.flags.hasRef
                                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit;
                            end
                            if MRSCont.flags.hasWater
                                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit;
                            end
                        else
                            names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk});
                             for f = 1 : length(names)
                                 if ~strcmp(names{f},'ppm')
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.(names{f});
                                 end
                             end
                              names = fields(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk});
                             for f = 1 : length(names)
                                 if ~strcmp(names{f},'ppm')
                                    MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.(names{f})= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.(names{f});
                                 end
                             end
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs/scale;
                            if MRSCont.flags.hasRef
                                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit;
                            end
                            if MRSCont.flags.hasWater
                                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                                MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit;
                            end
                        end
                end
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).C{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).C{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).D{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).D{1,kk}.specs/scale;
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.baseline/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.baseline/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.baseline/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.res/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.res/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.res/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff1{1,kk}.data/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_diff2{1,kk}.data/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).conc_sum{1,kk}.data/scale;
                        else
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.fit= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.fit/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.baseline/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.baseline/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.baseline= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.baseline/scale;
                           MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.data/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.data/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.data= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.data/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff1_diff1{1,kk}.res/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).diff2_diff2{1,kk}.res/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.res= MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).sum_sum{1,kk}.res/scale;
                        end
                        if MRSCont.flags.hasRef
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).ref_ref{1,kk}.fit/scale;
                        end
                        if MRSCont.flags.hasWater
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                            MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit =  MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).w_w{1,kk}.fit/scale;
                        end
                end

            else
                error('This script works only on fully processed data. Run the whole Osprey pipeline first. Seg/Coreg is not needed')
            end
        end
    end
else
    [MRSCont] = osp_fitInitialise(MRSCont);
    for rr = 1 : Voxels
        for kk = 1 : MRSCont.nDatasets
            if Voxels < 2
                scale                 = MRSCont.fit.scale{kk};
            else
                scale                 = MRSCont.fit.scale{kk};                        
            end
            msg = sprintf('Scaling data from dataset %d out of %d total datasets...', kk, MRSCont.nDatasets);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            fprintf([reverseStr, msg]);
            if MRSCont.flags.isGUI && isfield(progressText,'String')
                set(progressText,'String' ,sprintf('Scaling data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
                drawnow
            end
            if MRSCont.flags.isUnEdited 
                if MRSCont.flags.hasMM %re_mm
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).mm{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).mm{1,kk}.specs/scale; %re_mm
                end %re_mm
                if MRSCont.flags.hasRef
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                end
                if MRSCont.flags.hasWater
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                end
                MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
            end

            if MRSCont.flags.isMEGA
                    if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs/scale;
                        if MRSCont.flags.hasRef
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                        end
                        if MRSCont.flags.hasWater
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                        end
                    else                            
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).diff1{1,kk}.specs/scale;
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).sum{1,kk}.specs/scale;
                        if MRSCont.flags.hasRef
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                        end
                        if MRSCont.flags.hasWater
                            MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                        end
                    end
            end
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).A{1,kk}.specs/scale;
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).B{1,kk}.specs/scale;
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).C{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).C{1,kk}.specs/scale;
                    MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).D{1,kk}.specs= MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).D{1,kk}.specs/scale;
                    if MRSCont.flags.hasRef
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).ref{1,kk}.specs/scale;
                    end
                    if MRSCont.flags.hasWater
                        MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs =  MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).w{1,kk}.specs/scale;
                    end
            end
        end
    end
end
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
for i = 1 : MRSCont.nDatasets
    SepFileList{i} =  split(MRSCont.files{i}, filesep);
    subject{i} = [SepFileList{i}{end-1}]; % Create subject name list
end
if MRSCont.flags.hasStatfile % Has stat csv file
    statCSV = readtable(MRSCont.file_stat, 'Delimiter', ',','ReadVariableNames',1); % Load it
    name = statCSV.Properties.VariableNames;
    group_idx = find(strcmp(name,'group'));
    if isempty(group_idx) % No group supplied so create grand mean only
        MRSCont.overview.groups = ones(MRSCont.nDatasets,1);
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    else %Get grouping variable
        MRSCont.overview.groups = statCSV{:,group_idx}; 
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    end
    if  ~strcmp(name,'subject') % No subject names stored in the container
        
        if length(subject)>1 && ~strcmp(subject{1},subject{2}) % ADd names according to BIDS with subfolder
            if ~strcmp(name,'subject')
                statCSV.subject = subject';
            end
        else
            if ~strcmp(name,'subject') % Add whole path as BIDS wasn't set up properly
                statCSV.subject = MRSCont.files';
            end
        end
        writetable(statCSV,[MRSCont.outputFolder  filesep  'stat.csv']); % Write names into csv files
    end

else % No csv file supplied
    MRSCont.overview.groups = ones(MRSCont.nDatasets,1); %Create a single group
    MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    statCSV = array2table(MRSCont.overview.groups,'VariableNames',{'group'});
    if length(subject)>1 && ~strcmp(subject{1},subject{2}) %Add names to the csv file
        statCSV.subject = subject';
    else
        statCSV.subject = MRSCont.files';
    end
    writetable(statCSV,[MRSCont.outputFolder  filesep  'stat.csv']);
end

% Set up the different group names
MRSCont.overview.groupNames = cell(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    MRSCont.overview.groupNames{g} = ['Group ' num2str(g)];
end

%Exclude datasets based on the exclude field in the MRSConainer. THis can
%be triggered by pressing the left (remove) and right (add) arrow buttons
%in the listbox of the GUI
if isfield(MRSCont, 'exclude')
    if~isempty(MRSCont.exclude)
        MRSCont.overview.groups(MRSCont.exclude) = 0;
    end
end

% Sort the spectra according to the groups
for rr = 1 : Voxels
    for ss = 1 : NoSubSpec % Loop over subspectra
        for g = 1 : MRSCont.overview.NoGroups % loop over groups
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(['g_' num2str(g)]).(SubSpecNames{ss}) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
        end
        MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).GMean.(SubSpecNames{ss}) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss})(1,MRSCont.overview.groups > 0);
    end
end

if MRSCont.flags.didFit
    for rr = 1 : Voxels
        % Sort the models according to the groups
        for sf = 1 : NoFit % loop over fits
            for g = 1 : MRSCont.overview.NoGroups % loop over groups
                MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(['g_' num2str(g)]).([FitNames{sf} '_' dataPlotNames{sf}]) = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}])(1,MRSCont.overview.groups == g);
            end
            MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).GMean.([FitNames{sf} '_' dataPlotNames{sf}]) = MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).([FitNames{sf} '_' dataPlotNames{sf}])(1,MRSCont.overview.groups > 0);
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
        MRSCont.overview.corr.Meas{cor} = statCSV{:,cor};
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
    for ss = 1 : NoSubSpec %loop over subspectra
        names = fields(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]));
        for g = 1 : length(names) % loop over groups
            tempSubSpec = zeros(length(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(SubSpecNames{ss})),MRSCont.info.(SubSpecNames{ss}).max_ndatapoint);
            for kk = 1 : length(MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(SubSpecNames{ss})) % Loop over datasets to generate a matrix
                try
                    tempSubSpec(kk,:) = MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(SubSpecNames{ss}){1,kk}.specs;
                catch
                    tempSubSpec(kk,:) = ones(1,MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{1}){1,1}.sz(1)) *nan;
                end
            end
            %Calculate mean and SD
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(['mean_' SubSpecNames{ss}]) = nanmean(real(tempSubSpec),1);
            MRSCont.overview.Osprey.(['sort_data_voxel_' num2str(rr)]).(names{g}).(['sd_' SubSpecNames{ss}]) = nanstd(real(tempSubSpec),1);
        end
        %Store ppm
        MRSCont.overview.Osprey.(['ppm_data_' SubSpecNames{ss}]) = MRSCont.overview.Osprey.(['all_data_voxel_' num2str(rr)]).(SubSpecNames{ss}){MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind}.ppm;
    end
end

if MRSCont.flags.didFit
    %Do the same for the models 
    for rr = 1 : Voxels
        for sf = 1 : NoFit %loop over fits
            names = fields(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]));
            for g = 1 : length(names) %Loop over groups
                    tempSubSpec = zeros(length(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}])),length(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1}.ppm));
                    tempSubRes = tempSubSpec;
                    tempSubdata = tempSubSpec;
                    tempSubBaseline = tempSubSpec;
                    tempInidivMetab = [];
                    for kk = 1 : length(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}])) % Loop over datasets to generate a matrices
                      tempSubSpec(kk,:) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit; %Fits
                      tempSubRes(kk,:) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res; % Residuals
                      tempSubdata(kk,:) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data; % spectra
                      if ~(strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w')) %Is not water
                        tempSubBaseline(kk,:) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline; % Baseline
                        fits = fields(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}); % names of the basis functions
                         for f = 6 : length(fits) % loop over basis functions
                                tempInidivMetab.(fits{f})(kk,:)= MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(fits{f});
                         end
                      end
                    end
                    %Calculate mean and SD
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubSpec),1);
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['sd_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubSpec),1);
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_res_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubRes),1);
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['sd_res_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubRes),1);
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubdata),1);
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['sd_data_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubdata),1);

                    if ~(strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w')) %Is not water
                        MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_baseline_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubBaseline),1);
                        MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['sd_baseline_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubBaseline),1);
                        for f = 6 : length(fits) % loop over basis functions
                                MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_' fits{f} '_' FitNames{sf} '_' dataPlotNames{sf}]) = nanmean(real(tempInidivMetab.(fits{f})),1);
                                MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['sd_' fits{f} '_' FitNames{sf} '_' dataPlotNames{sf}]) = nanstd(real(tempInidivMetab.(fits{f})),1);
                        end
                    end

                    %Store ppm
                    MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,1}.ppm;
            end
        end
    end

    %Make sure the means are aligned
    for rr = 1 : Voxels
        for sf = 1 : NoFit %Loop over fits
            names = fields(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]));
            for g = 1 : length(names) %loop over groups
                    if MRSCont.flags.isUnEdited
                        if strcmp(FitNames{sf}, 'off')
                            %Find the ppm of the maximum peak magnitude within the given range:
                            ppmindex=find(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1)==max(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1);
                            ppmmax=ppmrange(ppmindex);
                            refShift=(ppmmax-2.013);
                        else
                            refShift = 0;
                        end
                    end
                    if MRSCont.flags.isMEGA
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            ppmindex=find(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1)==max(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1);
                            ppmmax=ppmrange(ppmindex);
                            refShift=(ppmmax-2.013);
                        else
                            naa = 1;
                            ppmindex=find(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_off_A(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A<2.1)==max(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_off_A(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_off_A<2.1);
                            if isempty(ppmindex)
                                ppmindex=find(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1)==max(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.data(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1)));
                                ppmrange=MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm>2.9 & MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]).off_A{1,kk}.ppm<3.1);
                                naa = 0;
                            end
                            ppmmax=ppmrange(ppmindex);
                            if length(ppmmax) > 1
                                refShift = 0;
                            else if naa == 1
                                    refShift=(ppmmax-2.013);
                                else
                                    refShift=(ppmmax-3.02);

                            end
                            end
                        end
                    end
                    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                        if isfield(MRSCont.overview.Osprey.(['all_models_voxel_' num2str(rr)]), 'conc_diff1')
                            ppmindex=find(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1)==max(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_conc_sum<2.1);
                            ppmmax=ppmrange(ppmindex);
                            refShift=(ppmmax-2.013);
                        else
                            ppmindex=find(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_sum_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum<2.1)==max(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).mean_data_sum_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum<2.1)));
                            ppmrange=MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum(MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).ppm_fit_sum_sum<2.1);
                            ppmmax=ppmrange(ppmindex);
                            refShift=(ppmmax-2.013);
                        end
                    end
                     MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) = MRSCont.overview.Osprey.(['sort_fit_voxel_' num2str(rr)]).(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) - refShift;
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

end
