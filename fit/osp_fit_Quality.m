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
    OrderNamesFit = horzcat(OrderNamesFit, 'w');
end

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
end


for ss = 1 :NoFitSpecNames %Loop over fitted spectra

for sf = 1 : size(FitSpecNamesStruct.(FitSpecNames{ss}),2) %Loop over all fits
    bf = 1;
        if ~isempty(FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf})
            for kk = 1 : MRSCont.nDatasets(1) %Loop over all datasets
                switch MRSCont.opts.fit.method %Which model was used
                case 'Osprey'
                    if ~strcmp(FitSpecNames{ss}, 'ref') && ~strcmp(FitSpecNames{ss}, 'w') && ~strcmp(FitSpecNames{ss}, 'mm') % metabolite only                        
                        fitRangePPM = MRSCont.opts.fit.range;

                        dataToPlot  = op_takesubspec(MRSCont.processed.(FitSpecNames{ss}){kk},find(strcmp(MRSCont.processed.(FitSpecNames{ss}){kk}.names,FitSpecNamesStruct.(FitSpecNames{ss}){bf,sf})));
                        basisSet    = MRSCont.fit.resBasisSet.(FitSpecNames{ss}).(['np_sw_' num2str(dataToPlot.sz(1)) '_' num2str(dataToPlot.spectralwidth)]){bf,sf};
                        fitParams   = MRSCont.fit.results.(FitSpecNames{ss}).fitParams{bf,kk,sf};
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
                        
                    end
                case 'LCModel'
                    dataToPlot  = MRSCont.processed.(FitSpecNames{ss}){kk};
                    fitParams   = MRSCont.fit.results.(FitSpecNames{ss}).fitParams{bf,kk,sf};

                    % Get the LCModel plots we previously extracted from .coord
                    % etc.
                    [ModelOutput] = fit_LCModelParamsToModel(fitParams);

                end
                
                %NOW FIND THE STANDARD DEVIATION OF THE NOISE:
                noisewindow=dataToPlot.specs(dataToPlot.ppm>-2 & dataToPlot.ppm<0)./MRSCont.fit.scale{kk};
                ppmwindow2=dataToPlot.ppm(dataToPlot.ppm>-2 & dataToPlot.ppm<0)';

                P=polyfit(ppmwindow2,noisewindow,2);
                noise=noisewindow-polyval(P,ppmwindow2); 

                MRSCont.QM.relAmpl.(FitSpecNamesStruct.(FitSpecNames{ss}){1,sf})(kk) = sum(ModelOutput.residual.^2)/(std(real(noise))^2 * length(ModelOutput.residual));
                    
            end
        end
end
end
        
    


