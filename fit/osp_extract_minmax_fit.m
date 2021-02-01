function MRSCont = osp_extract_minmax_fit(MRSCont, which_spec, conc)
%% out = osp_extract_minmax_fit(MRSCont, kk, which_spec, conc)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%                                 'mm' re_mm
%       conc      = flag if concatenate fitting was used
%       stagFlag  = flag to decide whether basis functions should be plotted
%                   vertically staggered or simply over one another
%                   (optional)
%                    OPTIONS:   1 = staggered (default)
%                               0 = not staggered
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided

if nargin<3
    conc = 'diff1'; 
    if nargin < 2
        which_spec = 'off';
        if nargin<1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
for kk = 1 : MRSCont.nDatasets
    if (MRSCont.flags.isPRIAM == 1)
        if  strcmp(which_spec, 'conc')
            dataToScale{1}=op_takeVoxel(MRSCont.processed.(conc){kk},1);
            dataToScale{2}=op_takeVoxel(MRSCont.processed.(conc){kk},2);
        else
            if strcmp(which_spec, 'off')
                dataToScale{1}=op_takeVoxel(MRSCont.processed.A{kk},1);
                dataToScale{2}=op_takeVoxel(MRSCont.processed.A{kk},2);
            else
                dataToScale{1}=op_takeVoxel(MRSCont.processed.(conc){kk},1);
                dataToScale{2}=op_takeVoxel(MRSCont.processed.(conc){kk},2);
            end
        end

        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet{1}.(which_spec).water{MRSCont.info.(which_spec).unique_ndatapoint_indsort(kk)};
        else if strcmp(which_spec, 'conc')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet{1}.(which_spec){MRSCont.info.diff1.unique_ndatapoint_indsort(kk)};
            else if strcmp(which_spec, 'off')
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{1}.(which_spec){kk};
                else
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{1}.(which_spec){kk};
                end
            end
        end
    else
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec).water{MRSCont.info.(which_spec).unique_ndatapoint_indsort(kk)};
        else if strcmp(which_spec, 'conc')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){MRSCont.info.diff1.unique_ndatapoint_indsort(kk)};
            else if strcmp(which_spec, 'off')
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){kk};
                else
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){kk};
                end
            end
        end        
    end




    % Get the fit parameters
    if (MRSCont.flags.isPRIAM == 1)
        fitParams{1}   = MRSCont.fit.results{1}.(which_spec).fitParams{kk};
        fitParams{2}   = MRSCont.fit.results{2}.(which_spec).fitParams{kk};
    else
        fitParams   = MRSCont.fit.results{VoxelIndex(1), VoxelIndex(2)}.(which_spec).fitParams{kk};
    end


    for rr = 1 : size(fitParams,2)
        % Pack up into structs to feed into the reconstruction functions
        inputData.dataToFit                 = dataToScale{rr};
        inputData.basisSet                  = basisSet;

        if (MRSCont.flags.isPRIAM == 1)
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
        inputSettings.concatenated.Subspec  = conc;

        %%% 3. PREPARE LINES TO DISPLAY %%%
        % Extract data, ppm axes, fit, residual, baseline, and individual
        % metabolite contributions.
        % Obviously, depending on the fit method that has been used, we need to
        % reconstruct the plots with fit method functions
        switch fitMethod
            % Depending on whether metabolite or water data are to be
            % displayed, create the plots via different models
            case 'Osprey'
                if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                    % if water, use the water model
                    [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams{rr});
                else
                    % if metabolites, use the metabolite model
                    if strcmp(inputSettings.fitStyle,'Concatenated')
                        [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams{rr});
                    else
                        [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams{rr});
                    end
                end
            case 'OspreyAsym'
                if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                    % if water, use the water model
                    [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams{rr});
                else
                    % if metabolites, use the metabolite model
                    if strcmp(inputSettings.fitStyle,'Concatenated')
                        [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams{rr});
                    else
                        [ModelOutput] = fit_OspreyAsymParamsToModel(inputData, inputSettings, fitParams{rr});
                    end
                end                
            case 'OspreyNoLS'
                if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
                    % if water, use the water model
                    [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams{rr});
                else
                    % if metabolites, use the metabolite model
                    if strcmp(inputSettings.fitStyle,'Concatenated')
                        [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams{rr});
                    else
                        [ModelOutput] = fit_OspreyNoLSParamsToModel(inputData, inputSettings, fitParams{rr});
                    end
                end        
        end
        dataToPlot  = ModelOutput.data;
        fit         = ModelOutput.completeFit;

        % Determine a positive stagger to offset data, fit, residual, and 
        % baseline from the individual metabolite contributions
        stagData(rr) = 0.1*(max(abs(min(dataToPlot)), abs(max(dataToPlot))));
        maxPlot(rr) = max(dataToPlot + abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData(rr);
    end

    [~,ind]= max(maxPlot);
    if  ~strcmp(which_spec, 'conc')
        MRSCont.plot.fit.(which_spec).maxPlot(kk) = maxPlot(ind);
        MRSCont.plot.fit.(which_spec).stagData(kk) = stagData(ind);
    else
        MRSCont.plot.fit.(conc).maxPlot(kk) = maxPlot(ind);
        MRSCont.plot.fit.(conc).stagData(kk) = stagData(ind);        
    end
end
end
   