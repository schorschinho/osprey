% osp_fit_Quality.m
%   This function calculates the relative amplitude of the residual
%   compared to the standard deviation of the noise. This is one of the
%   seven quality control parameters defined in the MRS consensus paper by
%   Wilson et al. ( https://doi.org/10.1002/mrm.27742). The original
%   desciption can be found in Barros & Slotboom
%   (https://doi.org/10.1016/j.ab.2017.01.017)
%
%   USAGE:
%       [MRSCont]=osp_fit_Quality(MRSCont);
%
%   INPUTS:
%       dataToFit     = data which has been fitted.
%       fitParams     = output of the model
%       basisSet      = basisSet 
%
%   OUTPUTS:
%       MRSCont     = MRS container
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
%       2020-06-26: First version of the code.

function [MRSCont]=osp_fit_Quality(MRSCont)

%%% 2. INITIALIZE VARIABLES %%%
%Getting the names of the SubSpectra and Fits
FitNames = fieldnames(MRSCont.fit.results);
NoFit = length(fieldnames(MRSCont.fit.results));
dataPlotNames = FitNames;
tempFitNames = FitNames;
shift = 0;

%Getting the final model names (needed for concatenated fits)
for sf = 1 : NoFit
    switch MRSCont.opts.fit.method
        case {'Osprey', 'LCModel'}
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
        case 'OspreyNoLS'
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
    end
end
FitNames = tempFitNames;
NoFit = length(FitNames);

% Apply the same steps to the fits
for sf = 1 : NoFit %Loop over all fits
    for kk = 1 : MRSCont.nDatasets %Loop over all datasets
        switch MRSCont.opts.fit.method %Which model was used
            
            case 'Osprey'
                if ~(strcmp((FitNames{sf}), 'ref') || strcmp((FitNames{sf}), 'w')) % Not the water model           
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}){kk};
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
                end
                
           case 'OspreyNoLS'
                if ~(strcmp((FitNames{sf}), 'ref') || strcmp((FitNames{sf}), 'w')) % Not the water model           
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}){kk};
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
                        [ModelOutput] = fit_OspreyNoLSParamsToModel(inputData, inputSettings, fitParams);
                    end            
                end
                
            case 'LCModel'
                dataToPlot  = MRSCont.processed.(dataPlotNames{sf}){kk};
                fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
                [ModelOutput] = fit_LCModelParamsToModel(fitParams);
                
        end
        
    %NOW FIND THE STANDARD DEVIATION OF THE NOISE:
    noisewindow=dataToPlot.specs(dataToPlot.ppm>-2 & dataToPlot.ppm<0)./MRSCont.fit.scale{kk};
    ppmwindow2=dataToPlot.ppm(dataToPlot.ppm>-2 & dataToPlot.ppm<0)';

    P=polyfit(ppmwindow2,noisewindow,2);
    noise=noisewindow-polyval(P,ppmwindow2); 
    
    MRSCont.QM.relAmpl.(dataPlotNames{sf})(kk) = sum(ModelOutput.residual.^2)/(std(real(noise))^2 * length(ModelOutput.residual));
    
    end
end


