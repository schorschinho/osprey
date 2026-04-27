function[out] = osp_plotModelvariation(MRSCont, Dataset, SubSpec, Ex)
%% function[out] = osp_plotModelvariation(MRSCont, Dataset, SubSpec, Ex)
%
% Description: Generates "hot spot" plots that visualize the regions of
%              model divergence by examining the standard deviation across
%              them for each PPM bin.
%
% Input:     MRSCont = Osprey container
%            Metab = specific metabolite on which to perform SCA (default = 'tNAA')
%            Dataset = The index of the spectrum (or spectra) to view
%            SubSpec = Sub spectrum
%            Ex = Experiment
%
% Example usage:
%                   osp_plotSpecificationcurve(MRSCont, 'GABA', 'tCr', 2);
%                   (plots the GABA/tCr ratio from the GABA-edited difference spectrum)
%                   
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
MRSCont = [];
Dataset = 0;
SubSpec {mustBeInteger}= 1;
Ex {mustBeInteger}= 1;
end
%% Handle errors and set up
if size(MRSCont.opts.fit.ModelProcedure.metab,1)==1
    error('You need to specify more than 1 model JSON in dim 1!')
elseif ~(MRSCont.flags.didFit)    
    error('Must run OspreyFit first!')
end

NModels = size(MRSCont.opts.fit.ModelProcedure.metab,1);
NPart = MRSCont.nDatasets(1);

if Dataset == 0
    Dataset = 1:NPart;
end

if MRSCont.flags.isGUI
    out = figure( 'Visible', 'off' );
else
    out = figure;
end

%% Retrieve the fits from MRSCont.fit.results

Shift = 0.15; % Shift between (normalized) plots
PPM = MRSCont.fit.results.metab{1,1,SubSpec,Ex,1}.Data.ppm;
FITS = zeros(size(MRSCont.fit.results.metab{1,1,SubSpec,Ex,1}.Model{end}.fit.fit,1),NModels); % Prealocate fits
Rng = MRSCont.fit.results.metab{1,1,SubSpec,Ex,1}.Options{end}.optimFreqFitRange;
for jj_models=1:NModels
    FITS_NPart = ones(size(MRSCont.fit.results.metab{1,1,SubSpec,Ex,1}.Model{end}.fit.fit,1),NPart) * nan; % Per-model FIT matrix
    Rng_jjmod = MRSCont.fit.results.metab{1,1,SubSpec,Ex,jj_models}.Options{end}.optimFreqFitRange;
    Rng = [min(Rng(1),Rng_jjmod(1)), max(Rng(2),Rng_jjmod(2))];
    for jj_part=1:length(Dataset)
        FITS_NPart(:,jj_part) = real(MRSCont.fit.results.metab{1,Dataset(jj_part),SubSpec,Ex,jj_models}.Model{end}.fit.fit);
    end
    FITS(:,jj_models) = median(FITS_NPart,2,'omitnan'); % Take median over datasets
end

%% Construct area coloring using standard deviation across models

FITS = FITS./max(abs(FITS(:))); % Normalize the fits

SD = std(FITS,[],2).';
SD(PPM<Rng(1) | PPM>Rng(2)) = 0;    % Remove SD outside the PPM range
SD = (SD - min(SD));                % Remove any systematic offsets

% Ascertain and fix YLims:
ymin = min(FITS(:,1)+Shift);
ymax = Shift*(NModels+1)+max(FITS(:));
ylim([ymin ymax])

if ~all(SD == 0 | isnan(SD)) % Only add hotspot overlay if there's something to show
    [X, Y] = meshgrid(PPM, [ymin ymax]);   % create 2-row mesh that matches ppm resolution and viewport
    C = repmat(SD, 2, 1);
    
    % Create surface plot
    bg = surface(X, Y, zeros(size(X)), C, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    % Specify the colorbar
    MapVals = 10;   % Number of bins with shading
    IgnoreVals = 5; % Number of bins without
    MinCol = 0.9;   % Minimum color of shaded areas (1=white)
    MaxCol = 0.5;   % Maximum color of shaded areas
    
    % Update color map
    colormap(gca,[repmat([1,1,1],IgnoreVals,1);[ones(MapVals,1),linspace(MinCol,MaxCol,MapVals)',linspace(MinCol,MaxCol,MapVals)']])
    alpha(bg, 'flat');
    bg.AlphaData = C;                % transparency matches SD
end

% Loop over and add model names alongside the plot
[~,ModelName] = fileparts(MRSCont.opts.fit.ModelProcedure.metab(:,SubSpec));
hold on
for jj_models=1:NModels
    plot(PPM,FITS(:,jj_models)+Shift*jj_models,'k')
    text(Rng(1),Shift*jj_models, ['  ',strrep(ModelName{jj_models},'_',' ')], 'HorizontalAlignment', 'left','FontSize',12);
    yline(jj_models*Shift)
end
set(gca, 'XDir','reverse','XLim',Rng, 'ytick',[], 'FontSize',12)

end