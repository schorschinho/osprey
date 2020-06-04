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
% Fall back to defaults if not provided
    if nargin<1
        error('ERROR: no input Osprey container specified.  Aborting!!');
    end

%%% 2. INITIALIZE VARIABLES %%%
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
FitNames = fieldnames(MRSCont.fit.results);
NoFit = length(fieldnames(MRSCont.fit.results));
dataPlotNames = FitNames;
tempFitNames = FitNames;
shift = 0;
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

%%% 3. INTERPOLATION & NORMALIZATION %%%
% Processed data
refProcessTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressbar = waitbar(0,'Start','Name','Osprey Overview');
    waitbar(0,progressbar,sprintf('Gathered spectra from subspectrum %d out of %d total subspectra...\n', 0, NoSubSpec))
end
MRSCont.overview.all_data = MRSCont.processed;
for ss = 1 : NoSubSpec
    msg = sprintf('Gathering spectra from subspectrum %d out of %d total subspectra...\n', ss, NoSubSpec);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));    
    for kk = 1 : MRSCont.nDatasets
        if MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1) < MRSCont.info.(SubSpecNames{ss}).max_ndatapoint
            ppmRangeData        = MRSCont.processed.(SubSpecNames{ss}){1,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind}.ppm';
            ppmRangeDataToInt       = MRSCont.processed.(SubSpecNames{ss}){1,kk}.ppm;
            ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
            if sum(ppmIsInDataRange) == 0
                ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
            end
            MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
            MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm = ppmRangeData;
            if mod(size(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),2)==0
                MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.fids=ifft(fftshift(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),[],MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t);
            else
                MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.fids=ifft(circshift(fftshift(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),1),[],MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t);
            end

        end
    end
    if MRSCont.flags.isGUI        
        waitbar(ss/NoSubSpec,progressbar,sprintf('Gathered spectra from subspectrum %d out of %d total subspectra...\n', ss, NoSubSpec))
    end
end


for ss = 1 : NoSubSpec
for kk = 1 : MRSCont.nDatasets
 %Find the ppm of the maximum peak magnitude within the given range:
    if MRSCont.flags.isUnEdited  
        %Find the ppm of the maximum peak magnitude within the given range:
        ppmindex=find(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm<2.1)==max(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm<2.1)));
        ppmrange=MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm>1.9 & MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm<2.1);
        ppmmax=ppmrange(ppmindex);
        refShift=(ppmmax-2.013);
    end
    if MRSCont.flags.isMEGA
        ppmindex=find(MRSCont.overview.all_data.sum{1,kk}.specs(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1)==max(MRSCont.overview.all_data.sum{1,kk}.specs(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1)));
        ppmrange=MRSCont.overview.all_data.sum{1,kk}.ppm(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1);
        ppmmax=ppmrange(ppmindex);
        refShift=(ppmmax-2.013);
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
        ppmindex=find(MRSCont.overview.all_data.sum{1,kk}.specs(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1)==max(MRSCont.overview.all_data.sum{1,kk}.specs(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1)));
        ppmrange=MRSCont.overview.all_data.sum{1,kk}.ppm(MRSCont.overview.all_data.sum{1,kk}.ppm>1.9 & MRSCont.overview.all_data.sum{1,kk}.ppm<2.1);
        ppmmax=ppmrange(ppmindex);
        refShift=(ppmmax-2.013);
    end    

 MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm = MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm - refShift;
end    
end

fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    pause(1)
    waitbar(0,progressbar,sprintf('Gathered fit models from fit %d out of %d total fits...\n', 0, NoFit))
end
reverseStr = '';
% Fits
for sf = 1 : NoFit
    msg = sprintf('Gathering fit models from fit %d out of %d total fits...\n', sf, NoFit);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
    for kk = 1 : MRSCont.nDatasets
        switch MRSCont.opts.fit.method
        case 'Osprey'
        if strcmp((FitNames{sf}), 'ref') || strcmp((FitNames{sf}), 'w')
            % if water, use the water model
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).water{MRSCont.info.(FitNames{sf}).unique_ndatapoint_indsort(kk)};
            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
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
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      = ModelOutput.ppm;
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;            
        else
            % if metabolites, use the metabolite model
            fitRangePPM = MRSCont.opts.fit.range;
            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}){MRSCont.info.A.unique_ndatapoint_indsort(kk)};
            dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
            % Get the fit parameters
            fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
            % Pack up into structs to feed into the reconstruction functions
            inputData.dataToFit                 = dataToPlot;
            inputData.basisSet                  = basisSet;
            inputSettings.scale                 = MRSCont.fit.scale{kk};
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
            if ~isnan(ModelOutput.completeFit)
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = ModelOutput.completeFit;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = ModelOutput.baseline;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  ModelOutput.ppm;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = ModelOutput.residual;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = ModelOutput.data;
                    
                    
                     if strcmp(FitNames{sf}, 'mm') %re_mm
                         for n = 1 : 4 + MRSCont.fit.basisSet.nMM
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' MRSCont.fit.basisSet.name{n}])  = ModelOutput.indivMets(:,n);                    
                        end
                        idx_NAA  = 4;
                        idx_Cr  = 1;
                        idx_CrCH2  = 2;
                        if ~isempty(idx_CrCH2)
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_CrCH2);
                        else
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) ;
                        end
                    if MRSCont.opts.fit.fitMM == 1
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,5:end),2);   
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitNAA  = sum(ModelOutput.indivMets(:,4),2);
                        
                    end
                    %section to write out MM_clean spectra
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.MM_clean = ModelOutput.data -sum(ModelOutput.indivMets(:,1:4),2);
                     else%re_mm
                        for n = 1 : MRSCont.fit.basisSet.nMets + MRSCont.fit.basisSet.nMM
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(['fit' MRSCont.fit.basisSet.name{n}])  = ModelOutput.indivMets(:,n);                    
                        end
                        idx_NAA  = find(strcmp(basisSet.name,'NAA'));
                        idx_NAAG  = find(strcmp(basisSet.name,'NAAG'));
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = ModelOutput.indivMets(:,idx_NAA) + ModelOutput.indivMets(:,idx_NAAG);
                        idx_Cr  = find(strcmp(basisSet.name,'Cr'));
                        idx_PCr  = find(strcmp(basisSet.name,'PCr'));
                        idx_CrCH2  = find(strcmp(basisSet.name,'CrCH2'));
                        if ~isempty(idx_CrCH2)
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_PCr)+ ModelOutput.indivMets(:,idx_CrCH2);
                        else
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = ModelOutput.indivMets(:,idx_Cr) + ModelOutput.indivMets(:,idx_PCr);
                        end
                        idx_1  = find(strcmp(basisSet.name,'GPC'));
                        idx_2  = find(strcmp(basisSet.name,'PCh'));
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCho  = ModelOutput.indivMets(:,idx_1) + ModelOutput.indivMets(:,idx_2);
                        idx_1  = find(strcmp(basisSet.name,'Glu'));
                        idx_2  = find(strcmp(basisSet.name,'Gln'));
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fitGlx  = ModelOutput.indivMets(:,idx_1) + ModelOutput.indivMets(:,idx_2);
                    if MRSCont.opts.fit.fitMM == 1
                        MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittMM  = sum(ModelOutput.indivMets(:,MRSCont.fit.basisSet.nMets+1:end),2);                    
                    end
                     end %re_mm
                else
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline      = nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm      =  nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res      = nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittNAA  = nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fittCr  = nan;
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = nan;
                end
        end
        end        
    end
    if MRSCont.flags.isGUI        
        waitbar(sf/NoFit,progressbar,sprintf('Gathered fit models from fit %d out of %d total fits...\n', sf, NoFit))
    end
end
fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    pause(1)
    waitbar(0,progressbar,sprintf('Interpolated fit models from fit %d out of %d total fits...\n', 0, NoFit))
end
reverseStr = '';

for sf = 1 : NoFit
    for kk = 1 : MRSCont.nDatasets
        temp_fit_sz.([FitNames{sf} '_' dataPlotNames{sf}])(1,kk)= length(MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit);
    end
    [max_point_fit.([FitNames{sf} '_' dataPlotNames{sf}]),max_ind_fit.([FitNames{sf} '_' dataPlotNames{sf}])] = max(temp_fit_sz.([FitNames{sf} '_' dataPlotNames{sf}]));
end


for sf = 1 : NoFit
    msg = sprintf('Interpolating fit models from fit %d out of %d total fits...\n', sf, NoFit);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));     
    for kk = 1 : MRSCont.nDatasets
        if length(MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit) < max_point_fit.([FitNames{sf} '_' dataPlotNames{sf}])
                    ppmRangeData        = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,max_ind_fit.([FitNames{sf} '_' dataPlotNames{sf}])}.ppm';
                    ppmRangeDataToInt       = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm;
                    ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                    if sum(ppmIsInDataRange) == 0
                        ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                    end
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');                   
                    if ~strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w')
                         MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');                        
                         names = fields(MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk});
                         for f = 6 : length(names)
                            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(names{f})= interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(names{f})(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                         end                                                
                    end
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm = ppmRangeData';
                    MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data-MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit;
        end
    end
    if MRSCont.flags.isGUI        
        waitbar(sf/NoFit,progressbar,sprintf('Interpolated fit models from fit %d out of %d total fits...\n', sf, NoFit))
    end
end

for sf = 1 : NoFit
    if ~strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w')
        for kk = 1 : MRSCont.nDatasets        
             %Find the ppm of the maximum peak magnitude within the given range:
             if MRSCont.flags.isUnEdited
                ppmindex=find(MRSCont.overview.all_models.off_A{1,kk}.data(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1)==max(MRSCont.overview.all_models.off_A{1,kk}.data(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1)));
                ppmrange=MRSCont.overview.all_models.off_A{1,kk}.ppm(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1);
             end
             if MRSCont.flags.isMEGA
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    ppmindex=find(MRSCont.overview.all_models.conc_sum{1,kk}.data(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.all_models.conc_sum{1,kk}.data(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1)));
                    ppmrange=MRSCont.overview.all_models.conc_sum{1,kk}.ppm(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1);
                else
                    ppmindex=find(MRSCont.overview.all_models.off_A{1,kk}.data(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1)==max(MRSCont.overview.all_models.off_A{1,kk}.data(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1)));
                    ppmrange=MRSCont.overview.all_models.off_A{1,kk}.ppm(MRSCont.overview.all_models.off_A{1,kk}.ppm>1.9 & MRSCont.overview.all_models.off_A{1,kk}.ppm<2.1);
                end                
             end
             if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    ppmindex=find(MRSCont.overview.all_models.conc_sum{1,kk}.data(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.all_models.conc_sum{1,kk}.data(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1)));
                    ppmrange=MRSCont.overview.all_models.conc_sum{1,kk}.ppm(MRSCont.overview.all_models.conc_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.conc_sum{1,kk}.ppm<2.1);
                else
                    ppmindex=find(MRSCont.overview.all_models.sum_sum{1,kk}.data(MRSCont.overview.all_models.sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.sum_sum{1,kk}.ppm<2.1)==max(MRSCont.overview.all_models.sum_sum{1,kk}.data(MRSCont.overview.all_models.sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.sum_sum{1,kk}.ppm<2.1)));
                    ppmrange=MRSCont.overview.all_models.sum_sum{1,kk}.ppm(MRSCont.overview.all_models.sum_sum{1,kk}.ppm>1.9 & MRSCont.overview.all_models.sum_sum{1,kk}.ppm<2.1);
                end           
             end
            ppmmax=ppmrange(ppmindex);
            MRSCont.overview.refShift(kk)=(ppmmax-2.013);
            MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.ppm - MRSCont.overview.refShift(kk);
        end
    end
end

fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    pause(1)
    waitbar(0,progressbar,sprintf('Scaled data from dataset %d out of %d total datasetss...\n', 0, MRSCont.nDatasets))
end
reverseStr = '';
%%% 3. SCALING DATA  %%%
for kk = 1 : MRSCont.nDatasets
    msg = sprintf('Scaling data from dataset %d out of %d total datasetss...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));       
    if isfield(MRSCont, 'quantify')
        if MRSCont.flags.isUnEdited            
            if MRSCont.flags.hasMM %re_mm
                MRSCont.overview.all_data.mm{1,kk}.specs =  MRSCont.overview.all_data.mm{1,kk}.specs/MRSCont.fit.scale{kk}; %re_mm
                MRSCont.overview.all_models.mm_mm{1,kk}.fit =  MRSCont.overview.all_models.mm_mm{1,kk}.fit/MRSCont.fit.scale{kk};
                MRSCont.overview.all_models.mm_mm{1,kk}.MM_clean =  MRSCont.overview.all_models.mm_mm{1,kk}.MM_clean/MRSCont.fit.scale{kk};
                MRSCont.overview.all_models.mm_mm{1,kk}.fitNAA =  MRSCont.overview.all_models.mm_mm{1,kk}.fitNAA/MRSCont.fit.scale{kk};
                MRSCont.overview.all_models.mm_mm{1,kk}.fittCr =  MRSCont.overview.all_models.mm_mm{1,kk}.fittCr/MRSCont.fit.scale{kk};
            end %re_mm
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.scale{kk};
                MRSCont.overview.all_models.ref_ref{1,kk}.fit =  MRSCont.overview.all_models.ref_ref{1,kk}.fit/MRSCont.fit.scale{kk};
            end
            if MRSCont.flags.hasWater
                MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.scale{kk};
                MRSCont.overview.all_models.w_w{1,kk}.fit =  MRSCont.overview.all_models.w_w{1,kk}.fit/MRSCont.fit.scale{kk};
            end
            names = fields(MRSCont.overview.all_models.off_A{1,kk});
             for f = 1 : length(names)
                 if ~strcmp(names{f},'ppm')
                    MRSCont.overview.all_models.off_A{1,kk}.(names{f})= MRSCont.overview.all_models.off_A{1,kk}.(names{f})/(MRSCont.fit.results.off.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.off.fitParams{1,kk}.ampl(idx_Cr));
                 end
             end
            MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.scale{kk};
        end

        if MRSCont.flags.isMEGA
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    MRSCont.overview.all_models.conc_diff1{1,kk}.fit= MRSCont.overview.all_models.conc_diff1{1,kk}.fit/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_sum{1,kk}.fit= MRSCont.overview.all_models.conc_sum{1,kk}.fit/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_diff1{1,kk}.baseline= MRSCont.overview.all_models.conc_diff1{1,kk}.baseline/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_sum{1,kk}.baseline= MRSCont.overview.all_models.conc_sum{1,kk}.baseline/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_diff1{1,kk}.res= MRSCont.overview.all_models.conc_diff1{1,kk}.res/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_sum{1,kk}.res= MRSCont.overview.all_models.conc_sum{1,kk}.res/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_diff1{1,kk}.data= MRSCont.overview.all_models.conc_diff1{1,kk}.data/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_models.conc_sum{1,kk}.data= MRSCont.overview.all_models.conc_sum{1,kk}.data/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.diff1{1,kk}.specs= MRSCont.overview.all_data.diff1{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.sum{1,kk}.specs= MRSCont.overview.all_data.sum{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};     
                    if MRSCont.flags.hasRef
                        MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                        MRSCont.overview.all_models.ref_ref{1,kk}.fit =  MRSCont.overview.all_models.ref_ref{1,kk}.fit/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    end
                    if MRSCont.flags.hasWater
                        MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                        MRSCont.overview.all_models.w_w{1,kk}.fit =  MRSCont.overview.all_models.w_w{1,kk}.fit/(MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_PCr)+ MRSCont.fit.results.conc.fitParams{1,kk}.ampl(idx_Cr));
                    end
                else
                    MRSCont.overview.all_models.off_A{1,kk}.fit= MRSCont.overview.all_models.off_A{1,kk}.fit/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.fit= MRSCont.overview.all_models.diff1_diff1{1,kk}.fit/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.off_A{1,kk}.baseline= MRSCont.overview.all_models.off_A{1,kk}.baseline/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.baseline= MRSCont.overview.all_models.diff1_diff1{1,kk}.baseline/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.off_A{1,kk}.res= MRSCont.overview.all_models.off_A{1,kk}.res/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.res= MRSCont.overview.all_models.diff1_diff1{1,kk}.res/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.off_A{1,kk}.data= MRSCont.overview.all_models.off_A{1,kk}.data/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.data= MRSCont.overview.all_models.diff1_diff1{1,kk}.data/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_Cr))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.diff1{1,kk}.specs= MRSCont.overview.all_data.diff1{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG))/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_data.sum{1,kk}.specs= MRSCont.overview.all_data.sum{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG))/MRSCont.fit.scale{kk};
                    if MRSCont.flags.hasRef
                        MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG))/MRSCont.fit.scale{kk};
                        MRSCont.overview.all_models.ref_ref{1,kk}.fit =  MRSCont.overview.all_models.ref_ref{1,kk}.fit/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    end
                    if MRSCont.flags.hasWater
                        MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG))/MRSCont.fit.scale{kk};
                        MRSCont.overview.all_models.w_w{1,kk}.fit =  MRSCont.overview.all_models.w_w{1,kk}.fit/(MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAA)+ MRSCont.fit.results.diff1.fitParams{1,kk}.ampl(idx_NAAG));
                    end
                end            
        end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.scale{kk};
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.scale{kk};
                MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.scale{kk};
                MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.scale{kk};
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    MRSCont.overview.all_models.conc_diff1{1,kk}.fit= MRSCont.overview.all_models.conc_diff1{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff2{1,kk}.fit= MRSCont.overview.all_models.conc_diff2{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_sum{1,kk}.fit= MRSCont.overview.all_models.conc_sum{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff1{1,kk}.baseline= MRSCont.overview.all_models.conc_diff1{1,kk}.baseline/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff2{1,kk}.baseline= MRSCont.overview.all_models.conc_diff2{1,kk}.baseline/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_sum{1,kk}.baseline= MRSCont.overview.all_models.conc_sum{1,kk}.baseline/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff1{1,kk}.res= MRSCont.overview.all_models.conc_diff1{1,kk}.res/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff2{1,kk}.res= MRSCont.overview.all_models.conc_diff2{1,kk}.res/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_sum{1,kk}.res= MRSCont.overview.all_models.conc_sum{1,kk}.res/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff1{1,kk}.data= MRSCont.overview.all_models.conc_diff1{1,kk}.data/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_diff2{1,kk}.data= MRSCont.overview.all_models.conc_diff2{1,kk}.data/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.conc_sum{1,kk}.data= MRSCont.overview.all_models.conc_sum{1,kk}.data/MRSCont.fit.scale{kk};                    
                else
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.fit= MRSCont.overview.all_models.diff1_diff1{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff2_diff2{1,kk}.fit= MRSCont.overview.all_models.diff2_diff2{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.sum_sum{1,kk}.fit= MRSCont.overview.all_models.sum_sum{1,kk}.fit/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.baseline= MRSCont.overview.all_models.diff1_diff1{1,kk}.baseline/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff2_diff2{1,kk}.baseline= MRSCont.overview.all_models.diff2_diff2{1,kk}.baseline/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.sum_sum{1,kk}.baseline= MRSCont.overview.all_models.sum_sum{1,kk}.baseline/MRSCont.fit.scale{kk}; 
                   MRSCont.overview.all_models.diff1_diff1{1,kk}.data= MRSCont.overview.all_models.diff1_diff1{1,kk}.data/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff2_diff2{1,kk}.data= MRSCont.overview.all_models.diff2_diff2{1,kk}.data/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.sum_sum{1,kk}.data= MRSCont.overview.all_models.sum_sum{1,kk}.data/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff1_diff1{1,kk}.res= MRSCont.overview.all_models.diff1_diff1{1,kk}.res/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.diff2_diff2{1,kk}.res= MRSCont.overview.all_models.diff2_diff2{1,kk}.res/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.sum_sum{1,kk}.res= MRSCont.overview.all_models.sum_sum{1,kk}.res/MRSCont.fit.scale{kk};
                end
                if MRSCont.flags.hasRef
                    MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.ref_ref{1,kk}.fit =  MRSCont.overview.all_models.ref_ref{1,kk}.fit/MRSCont.fit.scale{kk};
                end
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.scale{kk};
                    MRSCont.overview.all_models.w_w{1,kk}.fit =  MRSCont.overview.all_models.w_w{1,kk}.fit/MRSCont.fit.scale{kk};
                end            
        end
        
    else
        error('This script works only on fully processed data. Run the whole Osprey pipeline first. Seg/Coreg is not needed')
    end
    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Scaled data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets))
    end
end
fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    pause(1)
    close(progressbar);
end

%%% 4. SORTING DATA  %%%
if MRSCont.flags.hasStatfile
    statCSV = readtable(MRSCont.file_stat, 'Delimiter', ',','ReadVariableNames',1);
    name = statCSV.Properties.VariableNames;
    group_idx = find(strcmp(name,'group'));
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
        MRSCont.overview.sort_data.(['g_' num2str(g)]).(SubSpecNames{ss}) = MRSCont.overview.all_data.(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
    end
    MRSCont.overview.sort_data.GMean.(SubSpecNames{ss}) = MRSCont.overview.all_data.(SubSpecNames{ss})(1,:);
end

for sf = 1 : NoFit
    for g = 1 : MRSCont.overview.NoGroups
        MRSCont.overview.sort_fit.(['g_' num2str(g)]).([FitNames{sf} '_' dataPlotNames{sf}]) = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}])(1,MRSCont.overview.groups == g);
    end
    MRSCont.overview.sort_fit.GMean.([FitNames{sf} '_' dataPlotNames{sf}]) = MRSCont.overview.all_models.([FitNames{sf} '_' dataPlotNames{sf}])(1,:);
end





%%% 5. CALCULATE MEAN AND SD SPECTRA FOR VISUALIZATION %%%
%Data
for ss = 1 : NoSubSpec
    names = fields(MRSCont.overview.sort_data);
    for g = 1 : length(names)
        tempSubSpec = zeros(length(MRSCont.overview.sort_data.(names{g}).(SubSpecNames{ss})),MRSCont.overview.all_data.(SubSpecNames{1}){1,1}.sz(1));
        for kk = 1 : length(MRSCont.overview.sort_data.(names{g}).(SubSpecNames{ss}))
          tempSubSpec(kk,:) = MRSCont.overview.sort_data.(names{g}).(SubSpecNames{ss}){1,kk}.specs;
        end
        MRSCont.overview.sort_data.(names{g}).(['mean_' SubSpecNames{ss}]) = nanmean(real(tempSubSpec),1);
        MRSCont.overview.sort_data.(names{g}).(['sd_' SubSpecNames{ss}]) = nanstd(real(tempSubSpec),1);
    end
MRSCont.overview.(['ppm_data_' SubSpecNames{ss}]) = MRSCont.overview.all_data.(SubSpecNames{ss}){1,1}.ppm;
end    



%Fits
for sf = 1 : NoFit
    names = fields(MRSCont.overview.sort_fit);
    for g = 1 : length(names)
            tempSubSpec = zeros(length(MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1}),length(MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1}.ppm));
            tempSubBaseline = tempSubSpec;
            tempSubRes = tempSubSpec;
            tempSubdata = tempSubSpec;
            for kk = 1 : length(MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]))
              tempSubSpec(kk,:) = MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.fit;
              tempSubRes(kk,:) = MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.res;
              tempSubdata(kk,:) = MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.data;              
              if ~(strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w'))
                tempSubBaseline(kk,:) = MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.baseline;
                fits = fields(MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk});
                 for f = 6 : length(fits)
                        tempInidivMetab.(fits{f})(kk,:)= MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,kk}.(fits{f});
                 end
              end
            end
            MRSCont.overview.sort_fit.(names{g}).(['mean_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubSpec),1);
            MRSCont.overview.sort_fit.(names{g}).(['sd_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubSpec),1);
            MRSCont.overview.sort_fit.(names{g}).(['mean_res_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubRes),1);
            MRSCont.overview.sort_fit.(names{g}).(['sd_res_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubRes),1);
            MRSCont.overview.sort_fit.(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubdata),1);
            MRSCont.overview.sort_fit.(names{g}).(['sd_data_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubdata),1);            
            if ~(strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'ref_ref') || strcmp([FitNames{sf} '_' dataPlotNames{sf}], 'w_w'))
                MRSCont.overview.sort_fit.(names{g}).(['mean_baseline_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanmean(real(tempSubBaseline),1);
                MRSCont.overview.sort_fit.(names{g}).(['sd_baseline_' [FitNames{sf} '_' dataPlotNames{sf}]]) = nanstd(real(tempSubBaseline),1);
                for f = 6 : length(fits)
                        MRSCont.overview.sort_fit.(names{g}).(['mean_' fits{f} '_' FitNames{sf} '_' dataPlotNames{sf}]) = nanmean(real(tempInidivMetab.(fits{f})),1);
                        MRSCont.overview.sort_fit.(names{g}).(['sd_' fits{f} '_' FitNames{sf} '_' dataPlotNames{sf}]) = nanstd(real(tempInidivMetab.(fits{f})),1);
                end
            end
            MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) = MRSCont.overview.sort_fit.(names{g}).([FitNames{sf} '_' dataPlotNames{sf}]){1,1}.ppm;
    end
end
for sf = 1 : NoFit
    names = fields(MRSCont.overview.sort_fit);
    for g = 1 : length(names)
            if MRSCont.flags.isUnEdited  
                %Find the ppm of the maximum peak magnitude within the given range:
                ppmindex=find(MRSCont.overview.sort_fit.(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1)==max(MRSCont.overview.sort_fit.(names{g}).(['mean_data_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1)));
                ppmrange=MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])(MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])>1.9 & MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]])<2.1);
                ppmmax=ppmrange(ppmindex);
                refShift=(ppmmax-2.013);
            end
            if MRSCont.flags.isMEGA
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    ppmindex=find(MRSCont.overview.sort_fit.(names{g}).mean_data_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1)==max(MRSCont.overview.sort_fit.(names{g}).mean_data_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1)));
                    ppmrange=MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1);
                    ppmmax=ppmrange(ppmindex);
                    refShift=(ppmmax-2.013);
                else
                    ppmindex=find(MRSCont.overview.sort_fit.(names{g}).mean_data_off_A(MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A<2.1)==max(MRSCont.overview.sort_fit.(names{g}).mean_data_off_A(MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A<2.1)));
                    ppmrange=MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A(MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_off_A<2.1);
                    ppmmax=ppmrange(ppmindex);
                    refShift=(ppmmax-2.013);                  
                end
            end
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
                if isfield(MRSCont.overview.all_models, 'conc_diff1')
                    ppmindex=find(MRSCont.overview.sort_fit.(names{g}).mean_data_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1)==max(MRSCont.overview.sort_fit.(names{g}).mean_data_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1)));
                    ppmrange=MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_conc_sum<2.1);
                    ppmmax=ppmrange(ppmindex);
                    refShift=(ppmmax-2.013);
                else
                    ppmindex=find(MRSCont.overview.sort_fit.(names{g}).mean_data_sum_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum<2.1)==max(MRSCont.overview.sort_fit.(names{g}).mean_data_sum_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum<2.1)));
                    ppmrange=MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum(MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum>1.9 & MRSCont.overview.sort_fit.(names{g}).ppm_fit_sum_sum<2.1);
                    ppmmax=ppmrange(ppmindex);
                    refShift=(ppmmax-2.013);                  
                end
            end            
             MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) = MRSCont.overview.sort_fit.(names{g}).(['ppm_fit_' [FitNames{sf} '_' dataPlotNames{sf}]]) - refShift;
    end
end

%%% 6. READ CORRELATION DATA INTO THE STRUCT %%%
if MRSCont.flags.hasStatfile
    for cor = 1 : length(name)
        MRSCont.overview.corr.Names{cor} = name{cor};
    end
    for cor = 1 : length(name)
        MRSCont.overview.corr.Meas{cor} = statCSV{:,cor};
    end
end

%%% 7. CLEAN UP AND SAVE %%%
% Set exit flags and version
MRSCont.flags.didOverview          = 1;
MRSCont.ver.Over             = '100 Overview';

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;

if ~MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
%     save(fullfile(outputFolder, outputFile), 'MRSCont');
    MRSCont.flags.isGUI = 1;
else
%    save(fullfile(outputFolder, outputFile), 'MRSCont');
end

end

