function out = osp_plotGlobalConcentration(MRSCont,metabolite,quantification)
%% out = osp_plotGlobalConcentration(MRSCont,metabolite,quantification)
%   Creates a figure showing the quick integration maps of an MRSI scan
%
%   USAGE:
%       out = osp_plotGlobalConcentration(MRSCont,metabolite,quantification)
%
%   INPUTS:
%       out     = MATLAB figure handle
%       metabolite   = Target metabolite ('tNAA')
%       quantification = Which quantification to plot
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.
%%
% Fall back to defaults if not provided
if nargin<3
quantification = 'TissCorrWaterScaled';
    if nargin<2
        metabolite = 'tNAA';
        if nargin<1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end

%% Validate prerequisites
if ~MRSCont.flags.didOverview
    error('Trying to plot global concentrations, but Overview has not been performed yet. Run OspreyMRSIOverview first.')
end

%% Setup figure
out = figure;
set(out, 'Color', [1 1 1]); 
tiledlayout(2,4,'TileSpacing','compact')

% Global concentration bar plots
nexttile
X = categorical({'WM','GM'});  
b = bar(X,[MRSCont.GlobalConc.(quantification).(metabolite).C(2); MRSCont.GlobalConc.(quantification).(metabolite).C(1);], 'LineWidth', 2);
b.FaceColor = 'flat';
colors = [0.5 0.5 0.5; 1 1 1];
b.CData = colors;
hold on;
x = b.XEndPoints;  
errorbar(x, [MRSCont.GlobalConc.(quantification).(metabolite).C(2); MRSCont.GlobalConc.(quantification).(metabolite).C(1);],...
    [MRSCont.GlobalConc.(quantification).(metabolite).SE_C_WM MRSCont.GlobalConc.(quantification).(metabolite).SE_C_GM], 'k.', 'LineWidth', 1.5, 'CapSize', 10);
box off
set(gca,'TickDir','out')
ylabel([quantification ' concentration (i.u.)'])
title(['Global Concentrations ' metabolite],'interpreter','none')

% Tissue fraction historgram
nexttile
h = histogram(MRSCont.GlobalConc.(quantification).(metabolite).WM,20); hold on
h.FaceColor = [1 1 1];
h.EdgeColor = 'black';
h.LineWidth = 2;
h = histogram(MRSCont.GlobalConc.(quantification).(metabolite).GM,20);
h.FaceColor = [0 0 0];
h.FaceAlpha = 0.5;
h.EdgeColor = 'none';
legend('WM','GM')
box off
set(gca,'TickDir','out')
xlabel('tissue fraction')
ylabel('# of voxels')
title('Voxel Fraction Distribution')

% Linear regression plot
nexttile
scatter(MRSCont.GlobalConc.(quantification).(metabolite).WM,MRSCont.GlobalConc.(quantification).(metabolite).Q,'k');
box off
set(gca,'TickDir','out')
xlabel('WM fraction')
ylabel([quantification ' concentration (i.u.)'])
title(['WM vs concentration ' metabolite],'interpreter','none')

% Residual histogram
nexttile
histogram(MRSCont.GlobalConc.(quantification).(metabolite).residuals)
box off
set(gca,'TickDir','out')
xlabel([quantification ' concentration (i.u.)'])
ylabel('# of voxels')
title('Model Residual Distribution')

% Raw metabolite concentration map
nexttile
metab_image = squeeze(MRSCont.quantify.(quantification).(metabolite)(:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices));
metab_image(MRSCont.GlobalConc.(quantification).(metabolite).valid_mask_image == 0) = NaN;
imagesc(flip(rot90(metab_image),2))   
colormap gray
axis image
box off
axis off;
clim([MRSCont.GlobalConc.(quantification).(metabolite).residual_min MRSCont.GlobalConc.(quantification).(metabolite).Q_max]) 
set(gca,'TickDir','out')
title([metabolite ' Map'],'interpreter','none')

% Valid mask after applying all the thresholds
nexttile
imagesc(flip(rot90(squeeze(MRSCont.GlobalConc.(quantification).(metabolite).valid_mask_image)),2))
colormap gray
axis image
box off
axis off;
set(gca,'TickDir','out')   
title('Mask for Linear Regression')

% Predicted metabolite concentration map
nexttile
imagesc(flip(rot90(squeeze(MRSCont.GlobalConc.(quantification).(metabolite).Q_predicted_image)),2))
colormap gray
axis image
box off
axis off;
clim([MRSCont.GlobalConc.(quantification).(metabolite).residual_min MRSCont.GlobalConc.(quantification).(metabolite).Q_max])
set(gca,'TickDir','out')
title([metabolite ' Prediction Map'],'interpreter','none')

% Residual map
nexttile
imagesc(flip(rot90(squeeze(MRSCont.GlobalConc.(quantification).(metabolite).residuals_image)),2))
colormap gray
axis image
box off
axis off;
clim([MRSCont.GlobalConc.(quantification).(metabolite).residual_min MRSCont.GlobalConc.(quantification).(metabolite).Q_max])
set(gca,'TickDir','out')
cbar = colorbar;
cbar.Label.String = [quantification ' concentration (i.u.)'];
title('Model Residual Map')

end

   