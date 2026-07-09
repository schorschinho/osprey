function [MRSCont] = OspreyMRSI(jobFile,overwrite,stopAfterProcess)
%% MRSI analysis wrapper in Osprey
% This is script accompanies the beta version of the MRSI analysis pipeline
% in Osprey which was presented at ISMRM 2025. 
% 
% It will be part of the upcoming Osprey release 3.0.0
% Please make sure to remove any older Osprey versions from your Matlab
% path. Add the full OspreyMRSIbeta folder to the path. For data
% visualization you need to install FSL-eyes with the mrs-plugin. Currently
% you will also have to copy the viridis colourmap into the fsl folder.
% Find the location, e.g., ../FSLeyes/lib/python3.13/site-packages/fsleyes/assets/colourmaps
% copy the viridis color map file from osprey/mrsi/viridis.cmap into the
% folder, and add 'virdis   Viridis' to the order.txt file.
%
% The folder contains three example datasets from Philips and Siemens.
% Philips SPAR/SDAT and data/list files are most supported. All other data
% needs to be converted to nifti-mrs format using spec2nii ().
%
% Below you can find example function calls for the MRSI analysis
%% Parse input
if nargin < 3
    stopAfterProcess = 0;
    if nargin < 2
        overwrite = '00';
        if nargin<1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end
%% Initialization with OspreyJob
%
% Here data information (path to different MRSI files) are parsed into the
% MRSCont master structure. You can add voxel shifts and flips to the nifti
% output if needed
%
% MRSCont.opts.MRSI.nii_shifts = [0.5 0.5 0.5];
% MRSCont.opts.MRSI.nii_flip.cc = 1;

% Do you want interactive feedback during the processing?
% MRSCont.opts.MRSI.interactive = 0;

MRSCont = OspreyJob(jobFile,0,overwrite); % The which call is only need to ensure functinalty on other systems. You can directly put the full path to your job


%% Load data with OspreyLoad
%
% Vendor native raw data is parsed into the MRSCont including a fews
% preparatory steps (coil-combination, motion-correction, etc.) and
% exported as nifti-mrs files.

if ~MRSCont.flags.didLoadData
    MRSCont = OspreyLoad(MRSCont);
end

% Now we generate the quick maps for the raw data
field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{1});
for ff = 1 : length(field_names)
    MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{1}.(field_names{ff});
end

[MRSCont] = create_quickMaps(MRSCont);


%% Perform coregistration with OspreyCoreg
%
% MRSI volume is coregistered to the anatomical scan. This also generates
% the a voxel index mask at the resolution of the anatomical scan

if ~MRSCont.flags.didCoreg
    MRSCont = OspreyCoreg(MRSCont);
    if MRSCont.opts.MRSI.interactive
        out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI',1,0);
        uiwait(out);
        out = osp_plotInteractiveMRSIDataview(MRSCont);
        uiwait(out);
    end
end
%% Perform segmentation with OspreySeg
%
% The anatomical scan is segmented using SPM12. Here we also create
% the relevant masks for automated brain masking and quantification
%
% You can also add a square mask to accelerate processing
% Top right (confirm) corner has index [min(x) max(y)]
% MRSCont.opts.MRSI.outerMask.x = [2 33];
% MRSCont.opts.MRSI.outerMask.y = [4 42];
% MRSCont.opts.MRSI.outerMask.z = [1 5];

% Define threshold for lipid and brain mask here
% threshLipid is compared against the fractional amount of soft tissue if
% fLip > threshLipid lipid mask is set to 1
% threshBrain is compared against the sum of fWM and fGM if fWM + fGM >
% threshBrain brain mask is set to 1

% MRSCont.opts.MRSI.threshLipid = 0.1; % default is 0.1
% MRSCont.opts.MRSI.threshBrain = 0.3; % default is 0.3

if isempty(MRSCont.opts.MRSI.outerMask.x)   
    % If no mask is defined we define it here
    if isempty(MRSCont.opts.MRSI.outerMask.x)
        MRSCont.opts.MRSI.outerMask.x = [1,MRSCont.raw{1}.nXvoxels];
    end
    if isempty(MRSCont.opts.MRSI.outerMask.y)
        MRSCont.opts.MRSI.outerMask.y = [1,MRSCont.raw{1}.nYvoxels];
    end
    if isempty(MRSCont.opts.MRSI.outerMask.z)
        MRSCont.opts.MRSI.outerMask.z = [1,MRSCont.raw{1}.nZvoxels];
    end
    
    % Start interactive MRSI masking session 
    if  MRSCont.opts.MRSI.outerMask.Interactive
        out = osp_plotInteractiveMRSIMask(MRSCont, 'T1w_rMRSI');
        uiwait(out); 
        load(fullfile(MRSCont.outputFolder,MRSCont.outputFile))
    end
end

if ~MRSCont.flags.isPhantom
    if ~MRSCont.flags.didSeg   
        MRSCont = OspreySeg(MRSCont);
        if MRSCont.opts.MRSI.interactive
            out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI',1,1,1,1);
            uiwait(out);
        end
    end
end

%% Preprocess data with OspreyProcess
%
% The data is preprocessed prior to modeling. This includes eddy current
% correction (optional), removal of residual water (HSVD or L2-basis), 
% removal of lipid (L2-basis or L2-mask), and frequency correction.

% MRSCont.polResidNAA = 1; % Polarity NAA

% Options for auto phasing
% MRSCont.opts.MRSI.phase.type = 'none';
% MRSCont.opts.MRSI.phase.type = 'Cr-Cho';
% MRSCont.opts.MRSI.phase.type = 'auto_phase';
% MRSCont.opts.MRSI.phase.limits = [1.7,2.2];

% Options for MRSI nuisance signal removal are:
% MRSCont.opts.MRSI.NuisanceRemoval.water.type = 'none';
% MRSCont.opts.MRSI.NuisanceRemoval.water.type = 'HSVD';
% MRSCont.opts.MRSI.NuisanceRemoval.water.comp = 32; % Max number of components removed in HSVD
% MRSCont.opts.MRSI.NuisanceRemoval.water.type = 'L2-basis'; %L2 removal with simulated water basis functions
% MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments.Components = 2000;
% MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments.lineWidthRange= [10 80];
% MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments.PPMRange = [4.2 5.1];
% MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments.beta = 10;
% MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments.plotBasis = false;

% MRSCont.opts.MRSI.NuisanceRemoval.lipid.type = 'none';
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.type = 'L2-basis'; %L2 removal with simulated lipid basis functions
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.Components = 2000;
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.lineWidthRange= [10 80];
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.PPMRange = [0.3 1.900];
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.beta = 100;
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.plotBasis = false;
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.type = 'L2-mask'; %L2 removal with lipids defined by lip mask, requires segmentation!
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.Components = 32;
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.lineWidthRange= [10 80];
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.PPMRange = [0.3 1.85];
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.beta = 0.01;
% MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments.plotBasis = false;
% Options for MRSI frequency alignment are:
% MRSCont.opts.MRSI.FreqAlign.type = 'none'; % No frequency alignment

% MRSCont.opts.MRSI.FreqAlign.type = 'CC'; % Cross-correlation alignment of frequencies defined below
% MRSCont.opts.MRSI.FreqAlign.frequencies = [2.01,3.03,3.22];
% MRSCont.opts.MRSI.FreqAlign.polarity = [1,1,1];
% MRSCont.opts.MRSI.FreqAlign.lim = [1.85,6];
% MRSCont.opts.MRSI.FreqAlign.realpart = 0;
% MRSCont.opts.MRSI.FreqAlign.zerofill = 0;

% MRSCont.opts.MRSI.FreqAlign.type = 'CCwithLipRemoval'; % Cross-correlation alignment of frequencies defined below after Wavelet filter of baseline if lipid/noise > thresh
% MRSCont.opts.MRSI.FreqAlign.thresh = 10;
% MRSCont.opts.MRSI.FreqAlign.frequencies = [2.01,3.03,3.22];
% MRSCont.opts.MRSI.FreqAlign.polarity = [1,1,1];
% MRSCont.opts.MRSI.FreqAlign.lim = [1.85,4];
% MRSCont.opts.MRSI.FreqAlign.realpart = 0;
% MRSCont.opts.MRSI.FreqAlign.zerofill = 0;

% Flags for maximum echo processing
% MRSCont.opts.MRSI.MaxEcho.separate = 1;
% MRSCont.opts.MRSI.MaxEcho.tstart = 75; %When did the ADC start?
% MRSCont.opts.MRSI.MaxEcho.AdditionalPhasing = 1;
% MRSCont.opts.MRSI.MaxEcho.AdditionalFreqAlign = 1;

% Cross-correlation alignment of frequencies defined below
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.type = 'CC';
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.frequencies = [3.03,3.22,3.9];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.polarity = [1,1,1];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.lim = [2.5,4.5];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.realpart = 1;
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.zerofill = 1;

% Cross-correlation alignment of frequencies defined below after Wavelet filter of baseline if lipid/noise > thresh
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.type = 'CCwithLipRemoval'; 
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.thresh = 10;
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.frequencies = [2.01,3.03,3.22,3.9];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.polarity = [1,1,1,1];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.lim = [1.85,4.2];
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.realpart = 0;
% MRSCont.opts.MRSI.MaxEcho.FreqAlign.zerofill = 1;

if ~MRSCont.flags.didProcess
    MRSCont = OspreyProcess(MRSCont);
    if MRSCont.opts.MRSI.interactive
        out = osp_plotInteractiveMRSIDataview(MRSCont);
        uiwait(out);
    end
end
% Now we generate the quick maps for the processed data
field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{2});
for ff = 1 : length(field_names)
    MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{2}.(field_names{ff});
end
[MRSCont] = create_quickMaps(MRSCont);

if stopAfterProcess
    outputFolder = MRSCont.outputFolder;
    outputFile      = MRSCont.outputFile;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    
    if MRSCont.flags.isGUI
        MRSCont.flags.isGUI = 0;
        save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
        MRSCont.flags.isGUI = 1;
    else
       save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    end
    return
end
%% Linear-combination modeling
%
% Here linear-combination modeling is used to estimate metabolite
% amplitudes. The algorithm used here is a novel generalized algorithm
% which offers a lot of flexiblity in the analysis. The settings can be
% directly changed using so-called model-procedure json files. There are
% example model-procedure files in model-procedures folder.

if ~MRSCont.flags.didFit
    if isempty(MRSCont.opts.MRSI.MRSImask)
        MRSCont.MRSImask = flip(flip(squeeze(MRSCont.seg.tissue.brain(1,:,:,:)),1),3);
    else
        MRSCont.MRSImask = MRSCont.opts.MRSI.MRSImask;
    end
    
    MRSCont = MRSI_gLCM_wrapper(MRSCont, MRSCont.opts.MRSI.LCM.ModelProcedureFileMetabolites, MRSCont.opts.MRSI.LCM.ModelProcedureFileWater, ...
                                MRSCont.opts.MRSI.LCM.BasisSetFile, MRSCont.opts.MRSI.LCM.MetabSpecName, MRSCont.opts.MRSI.LCM.parallelComputing,MRSCont.opts.MRSI.LCM.zerofill);

    if MRSCont.opts.MRSI.interactive
        out = osp_plotInteractiveMRSIDataview(MRSCont);
        uiwait(out);
    end
end
%% Quantify and export results
%
% Here the amplitude estimates are translated into meaningful units of
% concentration. Current outputs include raw amplitudes, tCr referenced
% ratios, Water-Scaled, CSF-corrected Water-Scaled, and tissues and
% relaxation corrected molal concentrations (Gasparovic et al.).

% You can change the water visibility used for the Water-Scaled and
% CSF-corrected Water Scaled here
% MRSCont.opts.MRSI.WaterVisibility = 0.65; % Assuming pure white matter

% Automated QC is applied using 4 thresholds in the order SNR, FWHM, CRLB, and percentile. 
% First all voxels in the region are discarded if the have a lower SNR. 
% Second, voxels larger then the FWHM threshold (Hz) are discarded. These two
% thresholds are global and applied to all metabolites. Per metabolite
% based CRLB thresholding is applied next. Finally, a percentile based
% outlier rejection is applied. For example, percentile 99 means that the
% top 1% of the data is removed. 
MRSCont.opts.MRSI.Quantify.QC.SNRThreshold = 3;
MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold = 14;
MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold = 25;
MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold = 97;

if ~MRSCont.flags.didQuantify
    [MRSCont] = OspreyQuantifyMRSI(MRSCont,'metab');
end
%% Run additional analysis with Osprey Overview MRSI
%
% Here secondary analysis of the MRSI scan is performed. This includes the
% creation of representative gray and white matter spectra, calculation of
% global concentrations in gray and white matter using linear regression,
% and atlas based analysis of the MRSI results.

% Options for the calculations of the representative spectra the slice to
% include and the masking parameter of fGM + fWM per voxel. Only MRSI
% voxels greater then the threshold will be included
% MRSCont.opts.MRSI.RepSpectra.SliceIndices = [3];
% MRSCont.opts.MRSI.RepSpectra.fGMpfWM = 0.8;

% Options for the global concentrations are the slice to
% include and the masking parameter of fGM + fWM per voxel. Only MRSI
% voxels greater then the threshold will be included. Which quantification
% to use and which metabolites to include.
% MRSCont.opts.MRSI.GlobalConc.SliceIndices = [3];
% MRSCont.opts.MRSI.GlobalConc.fGMpfWM = 0.8;
% MRSCont.opts.MRSI.GlobalConc.metabolites = {'tNAA_Acetyl','tCr_methyl','tCho_mehtyl'};
% MRSCont.opts.MRSI.GlobalConc.quantities = {'TissCorrWaterScaled'};

% You can also visualize the atlas analysis in an interactive plot by
% calling osp_plotInteractiveAtlasAnalysis.

% Options for the atlas based analysis include which atlas to use with
% current options being 'neuromorphometrics' and 'AAL'. Additionally you
% can define whether left and right regions should also be combined in the
% statistics, the fractional threshold for a atlas, e.g., only MRSI voxels
% with that reach this threshold will be counted towards a specific label
% in the atlas. Further, 4 thresholds are applied for automated quality
% control in the order SNR, FWHM, CRLB, and SD. First all voxels in the
% region are discarded if the have a lower SNR. Second, voxels larger then the
% FWHM threshold (Hz) are discarded. Third, voxels that do not meet the
% CRLB threshold are discarded. You can define one CRLB per metabolite of
% interest. Lastly, mean and standard deviation of each region and
% metabolite are calculated and voxels that do not meet the SD threshold
% are discarded prior to the final statistical analysis. You can also
% define a list of metabolites to analyze and which quantification to use,
% inlcuding CRLBs as an option. The final statistics are stored as tsv
% files in the derivatives folder including per label statisistics as well
% as a summary of all MRSI voxels per label for further use.

% You can also visualize the atlas analysis in an interactive plot by
% calling osp_plotInteractiveAtlasAnalysis. 

% MRSCont.opts.MRSI.atlas.name = 'neuromorphometrics';
% MRSCont.opts.MRSI.atlas.CombineLR = 1;
% MRSCont.opts.MRSI.atlas.AtlasThreshold = 0.33;
% MRSCont.opts.MRSI.atlas.SNRThreshold = 3;
% MRSCont.opts.MRSI.atlas.FWHMThreshold = 14;
% MRSCont.opts.MRSI.atlas.CRLBThreshold = [20,20,20];
% MRSCont.opts.MRSI.atlas.SDThreshold = 3;
% MRSCont.opts.MRSI.atlas.metabolites = {'tNAA_Acetyl','tCr_methyl','tCho_mehtyl'};
% MRSCont.opts.MRSI.atlas.quantities = {'TissCorrWaterScaled','CRLBs'};

if ~MRSCont.flags.didOverview
    [MRSCont] = OspreyMRSIOverview(MRSCont);
    if MRSCont.opts.MRSI.interactive
        out = osp_plotInteractiveAtlasAnalysis(MRSCont,1,MRSCont.opts.MRSI.atlas.quantities{1},MRSCont.opts.MRSI.atlas.metabolites,'T1w_rMRSI','Pallidum_L',MRSCont.opts.MRSI.atlas.AtlasThreshold,1);
        uiwait(out);
        if MRSCont.flags.isUnEdited
            out = osp_plotGlobalConcentration(MRSCont,MRSCont.opts.MRSI.GlobalConc.metabolites{end},MRSCont.opts.MRSI.GlobalConc.quantities{1});
            uiwait(out);
        end       
        if MRSCont.flags.isMEGA
            out = osp_plotGlobalConcentration(MRSCont,MRSCont.opts.MRSI.GlobalConc.metabolites{end},MRSCont.opts.MRSI.GlobalConc.quantities{1});
            uiwait(out);
        end        
    end
end

%% Generate HTML Summary Report

% Here the semi-interactive HTML report is generated. There are a few
% optional parameters that can be passed here. For example the locations of
% the example spectra to be plotted and which quantifications/metabolites
% to show in the report.

% MRSCont.opts.MRSI.report.VoxelIndices = [round(MRSCont.raw{1}.nXvoxels/2),round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);
%                                         round(MRSCont.raw{1}.nXvoxels/2)+1,round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);
%                                         round(MRSCont.raw{1}.nXvoxels/2)+2,round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);];

% You can plot all the different quantification options depending on your
% input they include tCr, rawWaterScaled, CSFWaterScaled,
% TissCorrWaterScaled
% MRSCont.opts.MRSI.report.quantifications = {'tCr','rawWaterScaled','CSFrawWaterScaled','TissCorrWaterScaled'};
% 
% %Pick the metabolite names to be reported
% MRSCont.opts.MRSI.report.metabolites = {'tNAA_Acetyl','tCr_methyl','tCho_mehtyl'}; % 'tNAA','tCr', 'tCho', 'Glx','mI'

[MRSCont] = OspreyMRSIHTMLReport(MRSCont,1);
end