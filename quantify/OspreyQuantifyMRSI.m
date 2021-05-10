function [MRSCont] = OspreyQuantifyMRSI(MRSCont)
%% [MRSCont] = OspreyQuantifyMRSI(MRSCont)
%   This function transforms the raw amplitude parameters determined during
%   OspreyFit into MRSI amplitude maps.
%
%   By default, OspreyQuantify will report tCr ratios for all metabolites.
%   These values will not undergo any further correction for tissue
%   content, or relaxation.
%
%   If some sort of unsuppressed data has been provided, OspreyQuantify will
%   calculate concentration estimates in institutional units. These values
%   will have varying degrees of corrections and assumptions.
%
%   USAGE:
%       MRSCont = OspreyQuantifyMRSI(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Hege Zollner (Johns Hopkins University, 2021-01-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       202-01-06: First version of the code.


outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));
% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    msg = 'Trying to quantify data, but data have not been modelled yet. Run OspreyFit first.';
    fprintf(msg);
    error(msg);    
end

% Version check and updating log file
MRSCont.ver.Quant             = '1.0.0 Quant';
fprintf(['Timestamp %s ' MRSCont.ver.Osp '  ' MRSCont.ver.Quant '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));


%%% 0. CHECK WHICH QUANTIFICATIONS CAN BE DONE %%%
% tCr ratios can always be calculated (unless the unlikely case that Cr is
% not in the basis set, a case we'll omit for now).
qtfyCr = 1;

% Check which types of metabolite data are available
if MRSCont.flags.isUnEdited
    getResults = {'off'};
elseif MRSCont.flags.isMEGA
    if strcmpi(MRSCont.opts.fit.style, 'Separate')
        getResults = {'diff1', 'off'};
    elseif strcmpi(MRSCont.opts.fit.style, 'Concatenated')
        getResults = {'conc'};
    end
elseif MRSCont.flags.isHERMES
    if strcmpi(MRSCont.opts.fit.style, 'Separate')
        getResults = {'diff1', 'diff2', 'sum'};
    elseif strcmpi(MRSCont.opts.fit.style, 'Concatenated')
        getResults = {'conc'};
    end
elseif MRSCont.flags.isHERCULES
    if strcmpi(MRSCont.opts.fit.style, 'Separate')
        getResults = {'diff1', 'diff2', 'sum'};
    elseif strcmpi(MRSCont.opts.fit.style, 'Concatenated')
        getResults = {'conc'};
    end
end

% Check which types of water data are available
if [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [1 1]
        % If both water reference and short-TE water data have been
        % provided, use the one with shorter echo time.
        qtfyH2O     = 1;
        getResults{end+1} = 'w';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [1 0]
        % If only one type of water data has been provided, use it.
        qtfyH2O     = 1;
        getResults{end+1} = 'ref';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [0 1]
        % If only one type of water data has been provided, use it.
        qtfyH2O     = 1;
        getResults{end+1} = 'w';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [0 0]
        % If no water ref has been provided, only tCr ratios can be
        % provided.
        qtfyH2O     = 0;
end

% Check whether segmentation has been run, and whether tissue parameters
% exist. In that case, we can do CSF correction, and full tissue
% correction.
if qtfyH2O == 1 && MRSCont.flags.didSeg && isfield(MRSCont.seg, 'tissue') 
    qtfyCSF     = 1;
    qtfyTiss    = 1;
else
    qtfyCSF     = 0;
    qtfyTiss    = 0;
end

% Check whether tissue correction is available and whether GABA-edited
% MEGA-PRESS has been run. In this case, we can apply the alpha correction
% (Harris et al, J Magn Reson Imaging 42:1431-40 (2015)).
if qtfyTiss == 1 && MRSCont.flags.isMEGA && (strcmp(MRSCont.opts.editTarget{1},'GABA'))
    qtfyAlpha   = 1;
else if qtfyTiss == 1 && (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) && (strcmp(MRSCont.opts.editTarget{1},'GABA') || strcmp(MRSCont.opts.editTarget{2},'GABA')) 
    qtfyAlpha   = 1;
     else
        qtfyAlpha   = 0;
    end
end

warning('off','all');

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'QuantifyResults');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

%Find dimensionality of MRSI acquisition
dim = size(MRSCont.fit.results);

% Add combinations of metabolites to the basisset
if qtfyH2O
    water = 1;
else
    water = 0;
end

for ll = 1:length(getResults) - water
    if length(dim) == 2        
        MRSCont.quantify.metabs.(getResults{ll}) = MRSCont.fit.resBasisSet.(getResults{ll}){1,1}.name;
    else
        MRSCont.quantify.metabs.(getResults{ll}) = MRSCont.fit.resBasisSet.(getResults{ll}){1,1}.name;
    end
end

if qtfyH2O
    if length(dim) == 2        
        MRSCont.quantify.metabs.(getResults{end}) = MRSCont.fit.resBasisSet.(getResults{end}).water{1,1}.name;
    else
        MRSCont.quantify.metabs.(getResults{end}) = MRSCont.fit.resBasisSet.(getResults{end}).water{1,1}.name;
    end
end

for kk = 1:MRSCont.nDatasets
    for ll = 1:length(getResults)
        for mm = 1:length(MRSCont.quantify.metabs.(getResults{ll}))
            if length(dim) == 2
                for x = 1 : dim(1)
                    for y = 1 : dim(2)
                        MRSCont.quantify.amplMets{kk}.(getResults{ll}).(MRSCont.quantify.metabs.(getResults{ll}){mm})(x,y) = MRSCont.fit.results{x,y}.(getResults{ll}).fitParams{kk}.ampl(mm);
                    end % Y voxel
                end % X voxel
            else
                for z = 1 : dim(3)
                    for x = 1 : dim(1)
                        for y = 1 : dim(2)
                            MRSCont.quantify.amplMets{kk}.(getResults{ll})(x,y,z) = MRSCont.fit.results{x,y,z}.(getResults{ll}).fitParams{kk}.ampl;
                        end % Y voxel
                    end % X voxel  
                end % slices
            end
        end
    end
end
% MRSCont = addMetabComb(MRSCont, getResults);


%% Loop over all datasets
QuantifyTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyQuant',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 

    %%% 1. GET BASIS SET AND FIT AMPLITUDES %%%
%     metsName = MRSCont.quantify.metabs; % just for the names
%     amplMets = MRSCont.quantify.amplMets{kk};


%     %%% 2. GET CREATINE RATIOS %%%
%     % We can always do this, but let's use the flag just to be safe:
%     if qtfyCr
% 
%         % Extract metabolite amplitudes from fit
%         tCrRatios = quantCr(metsName, amplMets, getResults);
% 
%         % Save back to Osprey data container
%         for ll = 1:length(getResults)
%             MRSCont.quantify.(getResults{ll}).tCr{kk}  = tCrRatios.(getResults{ll});
%         end
% 
%     end
% 
% 
%     %%% 3. GET WATER-SCALED, TISSUE-UNCORRECTED RATIOS %%%
%     if qtfyH2O
%         amplWater = MRSCont.fit.results.(waterType).fitParams{kk}.ampl;
%         % Get repetition times
%         metsTR  = MRSCont.processed.A{kk}.tr * 1e-3;
%         waterTR = MRSCont.processed.(waterType){kk}.tr * 1e-3;
%         % Get echo times
%         metsTE  = MRSCont.processed.A{kk}.te * 1e-3;
%         waterTE = MRSCont.processed.(waterType){kk}.te * 1e-3;
%         % Calculate water-scaled, but not tissue-corrected metabolite levels
%         rawWaterScaled = quantH2O(metsName, amplMets, amplWater, getResults, metsTR, waterTR, metsTE, waterTE);
% 
%         % Save back to Osprey data container
%         for ll = 1:length(getResults)
%             MRSCont.quantify.(getResults{ll}).rawWaterScaled{kk} = rawWaterScaled.(getResults{ll});
%         end
% 
%     end
% 
% 
%     %%% 4. GET CSF CORRECTION %%%
%     if qtfyCSF
% 
%         % Apply CSF correction
%         fCSF = MRSCont.seg.tissue.fCSF(kk);
%         CSFWaterScaled = quantCSF(rawWaterScaled, fCSF, getResults);
% 
%         % Save back to Osprey data container
%         for ll = 1:length(getResults)        
%             MRSCont.quantify.(getResults{ll}).CSFWaterScaled{kk} = CSFWaterScaled.(getResults{ll});
%         end
%     end
% 
% 
%     %%% 5. GET TISSUE CORRECTION %%%
%     if qtfyTiss
% 
%         % Apply tissue correction
%         fGM = MRSCont.seg.tissue.fGM(kk);
%         fWM = MRSCont.seg.tissue.fWM(kk);
%         TissCorrWaterScaled = quantTiss(metsName, amplMets, amplWater, getResults, metsTR, waterTR, metsTE, waterTE, fGM, fWM, fCSF);
% 
%         % Save back to Osprey data container
%         for ll = 1:length(getResults)        
%             MRSCont.quantify.(getResults{ll}).TissCorrWaterScaled{kk} = TissCorrWaterScaled.(getResults{ll});
%         end
%     end
% 
% 
%     %%% 6. GET ALPHA CORRECTION (THIS IS FOR GABA ONLY AT THIS POINT) %%%
%     if qtfyAlpha
% 
%         % For now, this is done for GABA only; however, the principle
%         % could be extended to other metabolites, as long as we have some
%         % form of prior knowledge about their concentrations in GM and WM.
% 
%         % Calculate mean WM/GM fractions
%         meanfGM = mean(MRSCont.seg.tissue.fGM); % average GM fraction across datasets
%         meanfWM = mean(MRSCont.seg.tissue.fWM); % average WM fraction across datasets
%         [AlphaCorrWaterScaled, AlphaCorrWaterScaledGroupNormed] = quantAlpha(metsName,amplMets, amplWater,getResults, metsTR, waterTR, metsTE, waterTE, fGM, fWM, fCSF, meanfGM, meanfWM,MRSCont.opts.fit.coMM3);
% 
%         % Save back to Osprey data container
%         MRSCont.quantify.(getResults{1}).AlphaCorrWaterScaled{kk} = AlphaCorrWaterScaled;
%         MRSCont.quantify.(getResults{1}).AlphaCorrWaterScaledGroupNormed{kk} = AlphaCorrWaterScaledGroupNormed;
%     end
end
time = toc(QuantifyTime);
if MRSCont.flags.isGUI && isfield(progressText,'String')      
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf('... done.\n Elapsed time %f seconds\n',time);
MRSCont.runtime.Quantify = time;
%% Create tables
% Set up readable tables for each quantification.
% [MRSCont] = osp_createTable(MRSCont,'amplMets', getResults);
% if qtfyCr
%     [MRSCont] = osp_createTable(MRSCont,'tCr', getResults);
% end
% if qtfyH2O
%     [MRSCont] = osp_createTable(MRSCont,'rawWaterScaled', getResults);
% end
% if qtfyCSF
%     [MRSCont] = osp_createTable(MRSCont,'CSFWaterScaled', getResults);
% end
% if qtfyTiss
%     [MRSCont] = osp_createTable(MRSCont,'TissCorrWaterScaled', getResults);
% end
% if qtfyAlpha
%     [MRSCont] = osp_createTable(MRSCont,'AlphaCorrWaterScaled', getResults(1));
%     [MRSCont] = osp_createTable(MRSCont,'AlphaCorrWaterScaledGroupNormed', getResults(1));
% end
%% Clean up and save
% Set exit flags
MRSCont.flags.didQuantify           = 1;
diary off
% Save the metabolite tables as CSV structure
% exportCSV (MRSCont,saveDestination, getResults);

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
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

end

%%

%%% Add combinations of metabolites %%%
function MRSCont = addMetabComb(MRSCont, getResults)
%% Loop over all datasets
for kk = 1:MRSCont.nDatasets
    % tNAA NAA+NAAG
    for ll = 1:length(getResults)
        idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'NAA'));
        idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'NAAG'));
        if  ~isempty(idx_1) && ~isempty(idx_2)
            idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tNAA'));
            if isempty(idx_3)
                MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'tNAA';
            end
            idx_tNAA = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tNAA'));
            tNAA = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
            MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_tNAA) = tNAA;

        end
    end
    % Glx Glu+Gln
    for ll = 1:length(getResults)
        idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Glu'));
        idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Gln'));   
        if  ~isempty(idx_1) && ~isempty(idx_2)
            idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Glx'));
            if isempty(idx_3)
                MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'Glx';
            end
            idx_Glx = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Glx'));
            Glx = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
            MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_Glx) = Glx;
        end
    end
    % tCho GPC+PCh
    for ll = 1:length(getResults)
        idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GPC'));
        idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'PCh'));
        if  ~isempty(idx_1) && ~isempty(idx_2)
            idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tCho'));
            if isempty(idx_3)
                MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'tCho';
            end
            idx_tCho = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tCho'));
            tCho = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
            MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_tCho) = tCho;
        end
    end
    % tCr Cr+PCr
    for ll = 1:length(getResults)
        idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Cr'));
        idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'PCr'));
        if  ~isempty(idx_1) && ~isempty(idx_2)
            idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tCr'));
            if isempty(idx_3)
                MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'tCr';
            end
            idx_tCr = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'tCr'));
            tCr = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
            MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_tCr) = tCr;
        end
    end
    %Glc+Tau
    for ll = 1:length(getResults)
        idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Glc'));
        idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'Tau'));
        if  ~isempty(idx_1) && ~isempty(idx_2)
            idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GlcTau'));
            if isempty(idx_3)
                MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'GlcTau';
            end
            idx_GlcTau = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GlcTau'));
            GlcTau = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
            MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GlcTau) = GlcTau;
        end
    end
    %GABA+coMM3
    if strcmp(MRSCont.opts.fit.coMM3, '1to1GABA') % fixed GABA coMM3 model
        for ll = 1:length(getResults)        
            idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABA'));
            if  ~isempty(idx_1)
                idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                if isempty(idx_3)
                    MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'GABAplus';
                end
                idx_GABAp = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1);
                MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp) = GABAp;
            end
        end        
    else if strcmp(MRSCont.opts.fit.coMM3, '3to2MM') % fixed MM09 coMM3 model
             for ll = 1:length(getResults)
                idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABA'));
                idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'MM09'));
                if  ~isempty(idx_1) && ~isempty(idx_2)
                    idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                    if isempty(idx_3)
                        MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'GABAplus';
                    end
                    idx_GABAp = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                    GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
                    MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp) = GABAp;
                end
            end              
        else % Models with a separate comMM3 function or without a co-edited MM function         
            for ll = 1:length(getResults)
                idx_1 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABA'));
                idx_2 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'MM3co'));
                if  ~isempty(idx_1) && ~isempty(idx_2)
                    idx_3 = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                    if isempty(idx_3)
                        MRSCont.quantify.metabs.(getResults{ll}){length(MRSCont.quantify.metabs.(getResults{ll}))+1} = 'GABAplus';
                    end
                    idx_GABAp = find(strcmp(MRSCont.quantify.metabs.(getResults{ll}),'GABAplus'));
                    GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2);
                    MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp) = GABAp;
                end
            end
        end
    end
end
end
%%% /Calculate ratios to totale creatine %%%
%%%%%%%%%%%% BELOW ARE THE QUANTIFICATION FUNCTIONS %%%%%%%%%%%%

%%% Calculate ratios to totale creatine %%%
function tCrRatios = quantCr(metsName, amplMets, getResults)
metsName = metsName.(getResults{1});
% Calculate tCr ratios
idx_Cr  = find(strcmp(metsName,'Cr'));
idx_PCr = find(strcmp(metsName,'PCr'));
for ll = 1:length(getResults)
    if isempty(idx_Cr) && isempty(idx_PCr)
        error('Error in OspreyQuantify: Creatine ratios cannot be calculated because neither Cr nor PCr are included in the basis set.')
    elseif isempty(idx_Cr) && ~isempty(idx_PCr)
        tCr.(getResults{ll}) = amplMets.(getResults{ll})(idx_PCr);
    elseif ~isempty(idx_Cr) && isempty(idx_PCr)
        tCr.(getResults{ll}) = amplMets.(getResults{ll})(idx_Cr);
    elseif ~isempty(idx_Cr) && ~isempty(idx_PCr)   
        tCr.(getResults{ll}) = amplMets.(getResults{ll})(idx_Cr) + amplMets.(getResults{ll})(idx_PCr);
    end
end
% If separate fit of sub-spectra has been performed, normalize to 'off' or
% 'sum'
if isfield(tCr, 'off')
    tCrNorm = tCr.off;
elseif isfield(tCr, 'sum')
    tCrNorm = tCr.sum;
elseif isfield(tCr, 'conc')
    tCrNorm = tCr.conc;
end
for ll = 1:length(getResults)
    tCrRatios.(getResults{ll}) = amplMets.(getResults{ll})./tCrNorm;
end

end
%%% /Calculate ratios to totale creatine %%%



%%% Calculate raw water-scaled estimates %%%
function rawWaterScaled = quantH2O(metsName, amplMets, amplWater, getResults, metsTR, waterTR, metsTE, waterTE)

% Define constants
PureWaterConc       = 55500;            % mmol/L
WaterVisibility     = 0.65;             % assuming pure white matter
metsTE              = metsTE * 1e-3;    % convert to s
waterTE             = waterTE * 1e-3;   % convert to s
metsTR              = metsTR * 1e-3;    % convert to s
waterTR             = waterTR * 1e-3;   % convert to s

% Look up relaxation times

% Water
% From Wansapura et al. 1999 (JMRI)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
T1_Water            = 1.100;            % average of WM and GM, Wansapura et al. 1999 (JMRI)
T2_Water            = 0.095;            % average of WM and GM, Wansapura et al. 1999 (JMRI)

% Metabolites
for ll = 1:length(getResults)
    for kk = 1:length(metsName.(getResults{ll}))
        [T1_Metab_GM(kk), T1_Metab_WM(kk), T2_Metab_GM(kk), T2_Metab_WM(kk)] = lookUpRelaxTimes(metsName.(getResults{1}){kk});
        % average across GM and WM
        T1_Metab(kk) = mean([T1_Metab_GM(kk) T1_Metab_WM(kk)]);
        T2_Metab(kk) = mean([T2_Metab_GM(kk) T2_Metab_WM(kk)]);
        T1_Factor(kk) = (1-exp(-waterTR./T1_Water)) ./ (1-exp(-metsTR./T1_Metab(kk)));
        T2_Factor(kk) = exp(-waterTE./T2_Water) ./ exp(-metsTE./T2_Metab(kk));

        % Calculate
        rawWaterScaled.(getResults{ll})(kk,1) = (amplMets.(getResults{ll})(kk) ./ amplWater) .* PureWaterConc ...
            .* WaterVisibility .* T1_Factor(kk) .* T2_Factor(kk);
    end
end
end
%%% /Calculate raw water-scaled estimates %%%



%%% Calculate CSF-corrected water-scaled estimates %%%
function CSFWaterScaled = quantCSF(rawWaterScaled, fCSF, getResults)

% Simply divide the raw water-scaled, but tissue-uncorrected values by the
% non-CSF fraction:
for ll = 1:length(getResults)
    CSFWaterScaled.(getResults{ll}) = rawWaterScaled.(getResults{ll}) ./ (1 - fCSF);
end

end
%%% /Calculate CSF-corrected water-scaled estimates %%%



%%% Calculate tissue-corrected water-scaled estimates %%%
function TissCorrWaterScaled = quantTiss(metsName, amplMets, amplWater, getResults, metsTR, waterTR, metsTE, waterTE, fGM, fWM, fCSF)
% This function calculates water-scaled, tissue-corrected metabolite
% estimates in molal units, according to Gasparovic et al, Magn Reson Med
% 55:1219-26 (2006).

% Define Constants
% Water relaxation
% From Lu et al. 2005 (JMRI)
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)
T1w_WM    = 0.832;
T2w_WM    = 0.0792;
T1w_GM    = 1.331;
T2w_GM    = 0.110;
T1w_CSF   = 3.817;
T2w_CSF   = 0.503;

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al. 2006 (MRM) uses relative densities, ref to
% Ernst et al. 1993 (JMR)
% fGM = 0.78
% fWM = 0.65
% fCSF = 0.97
% such that
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84
concW_GM    = 43.30*1e3;
concW_WM    = 36.08*1e3;
concW_CSF   = 53.84*1e3;
molal_concW = 55.51*1e3;

% Gasparovic et al. method
% Calculate molal fractions from volume fractions (equivalent to eqs. 5-7 in Gasparovic et al., 2006)
molal_fGM  = (fGM*concW_GM) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
molal_fWM  = (fWM*concW_WM) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
molal_fCSF = (fCSF*concW_CSF) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);

% Metabolites
for ll = 1:length(getResults)
    for kk = 1:length(metsName.(getResults{ll}))
        [T1_Metab_GM(kk), T1_Metab_WM(kk), T2_Metab_GM(kk), T2_Metab_WM(kk)] = lookUpRelaxTimes(metsName.(getResults{1}){kk});
        % average across GM and WM
        T1_Metab(kk) = mean([T1_Metab_GM(kk) T1_Metab_WM(kk)]);
        T2_Metab(kk) = mean([T2_Metab_GM(kk) T2_Metab_WM(kk)]);

        % Calculate water-scaled, tissue-corrected molal concentration
        % estimates
        TissCorrWaterScaled.(getResults{ll})(kk,1) = (amplMets.(getResults{ll})(kk) ./ amplWater) .* molal_concW ...
            .* (molal_fGM  * (1 - exp(-waterTR/T1w_GM)) * exp(-waterTE/T2w_GM) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk))) + ...
                molal_fWM  * (1 - exp(-waterTR/T1w_WM)) * exp(-waterTE/T2w_WM) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk))) + ...
                molal_fCSF * (1 - exp(-waterTR/T1w_CSF)) * exp(-waterTE/T2w_CSF) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk)))) / ...
                (1 - molal_fCSF);
    end
end
end
%%% /Calculate CSF-corrected water-scaled estimates %%%



%%% Calculate alpha-corrected water-scaled GABA estimates %%%
function [AlphaCorrWaterScaled, AlphaCorrWaterScaledGroupNormed] = quantAlpha(metsName, amplMets, amplWater,getResults, metsTR, waterTR, metsTE, waterTE, fGM, fWM, fCSF, meanfGM, meanfWM,coMM3)
% This function calculates water-scaled, alpha-corrected GABA
% estimates in molal units, according to Gasparovic et al, Magn Reson Med
% 55:1219-26 (2006).

% Define Constants
% Water relaxation
% From Lu et al. 2005 (JMRI)
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)
T1w_WM    = 0.832;
T2w_WM    = 0.0792;
T1w_GM    = 1.331;
T2w_GM    = 0.110;
T1w_CSF   = 3.817;
T2w_CSF   = 0.503;

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al. 2006 (MRM) uses relative densities, ref to
% Ernst et al. 1993 (JMR)
% fGM = 0.78
% fWM = 0.65
% fCSF = 0.97
% such that
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84
concW_GM    = 43.30*1e3;
concW_WM    = 36.08*1e3;
concW_CSF   = 53.84*1e3;

% Calculate alpha correction factor for GABA
cWM = 1; % concentration of GABA in pure WM
cGM = 2; % concentration of GABA in pure GM
alpha = cWM/cGM;
CorrFactor = (meanfGM + alpha*meanfWM) / ((fGM + alpha*fWM) * (meanfGM + meanfWM));

% GABA (Harris et al, J Magn Reson Imaging 42:1431-1440 (2015))
idx_GABA  = find(strcmp(metsName.(getResults{1}),'GABA'));
[T1_Metab_GM, T1_Metab_WM, T2_Metab_GM, T2_Metab_WM] = lookUpRelaxTimes(metsName.(getResults{1}){idx_GABA});
% average across GM and WM
T1_Metab = mean([T1_Metab_GM T1_Metab_WM]);
T2_Metab = mean([T2_Metab_GM T2_Metab_WM]);
ConcIU_TissCorr_Harris = (amplMets.(getResults{1})(idx_GABA) ./ amplWater) ...
        .* (fGM * concW_GM * (1 - exp(-waterTR/T1w_GM)) * exp(-waterTE/T2w_GM) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)) + ...
            fWM * concW_WM * (1 - exp(-waterTR/T1w_WM)) * exp(-waterTE/T2w_WM) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)) + ...
            fCSF * concW_CSF * (1 - exp(-waterTR/T1w_CSF)) * exp(-waterTE/T2w_CSF) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)));

AlphaCorrWaterScaled = ConcIU_TissCorr_Harris / (fGM + alpha*fWM);
AlphaCorrWaterScaledGroupNormed = ConcIU_TissCorr_Harris * CorrFactor;

if ~strcmp(coMM3, 'none')
    % GABA (Harris et al, J Magn Reson Imaging 42:1431-1440 (2015))
    idx_GABAp  = find(strcmp(metsName.(getResults{1}),'GABAplus'));
    [T1_Metab_GM, T1_Metab_WM, T2_Metab_GM, T2_Metab_WM] = lookUpRelaxTimes(metsName.(getResults{1}){idx_GABA});
    % average across GM and WM
    T1_Metab = mean([T1_Metab_GM T1_Metab_WM]);
    T2_Metab = mean([T2_Metab_GM T2_Metab_WM]);
    ConcIU_TissCorr_Harris = (amplMets.(getResults{1})(idx_GABAp) ./ amplWater) ...
            .* (fGM * concW_GM * (1 - exp(-waterTR/T1w_GM)) * exp(-waterTE/T2w_GM) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)) + ...
                fWM * concW_WM * (1 - exp(-waterTR/T1w_WM)) * exp(-waterTE/T2w_WM) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)) + ...
                fCSF * concW_CSF * (1 - exp(-waterTR/T1w_CSF)) * exp(-waterTE/T2w_CSF) / ((1 - exp(-metsTR/T1_Metab)) * exp(-metsTE/T2_Metab)));

    AlphaCorrWaterScaled(2,1) = ConcIU_TissCorr_Harris / (fGM + alpha*fWM);
    AlphaCorrWaterScaledGroupNormed(2,1) = ConcIU_TissCorr_Harris * CorrFactor;
end
end
%%% /Calculate alpha-corrected water-scaled GABA estimates %%%



%%% Lookup function for metabolite relaxation times %%%
function [T1_GM, T1_WM, T2_GM, T2_WM] = lookUpRelaxTimes(metName)

% Look up table below
% T1 values for NAA, Glu, Cr, Cho, Ins from Mlynarik et al, NMR Biomed
% 14:325-331 (2001)
% T1 for GABA from Puts et al, J Magn Reson Imaging 37:999-1003 (2013)
% T2 values from Wyss et al, Magn Reson Med 80:452-461 (2018)
% T2 values are averaged between OCC and pACC for GM; and PVWM for WM
relax.Asc   = [1340 1190 (125+105)/2 172];
relax.Asp   = [1340 1190 (111+90)/2 148];
relax.Cr    = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting
relax.GABA  = [1310 1310 (102+75)/2 (102+75)/2]; % No WM estimate available; take GM estimate; both in good accordance with 88 ms reported by Edden et al
relax.Glc   = [1340 1190 (117+88)/2 155]; % Glc1: [1310 1310 (128+90)/2 156];
relax.Gln   = [1340 1190 (122+99)/2 168];
relax.Glu   = [1270 1170 (135+122)/2 124];
relax.Gly   = [1340 1190 (102+81)/2 152];
relax.GPC   = [1300 1080 (274+222)/2 218]; % This is the Choline singlet (3.21 ppm, tcho2 in the paper); glycerol is tcho: [1310 1310 (257+213)/2 182]; % choline multiplet is tcho1: [1310 1310 (242+190)/2 178];
relax.GSH   = [1340 1190 (100+77)/2 145]; % This is the cysteine signal (GSH1 in the paper), glycine is GSH: [1310 1310 (99+72)/2 145]; % glutamate is GSH2: [1310 1310 (102+76)/2 165];
relax.Lac   = [1340 1190 (110+99)/2 159];
relax.Ins   = [1230 1010 (244+229)/2 161];
relax.NAA   = [1470 1350 (253+263)/2 343]; % This is the 2.008 ppm acetyl signal (naa in the paper); aspartyl is naa1: [1310 1310 (223+229)/2 310];
relax.NAAG  = [1340 1190 (128+107)/2 185]; % This is the 2.042 ppm acetyl signal (naag in the paper); aspartyl is naag1: [1310 1310 (108+87)/2 180]; % glutamate is NAAG2: [1310 1310 (110+78)/2 157];
relax.PCh   = [1300 1080 (274+221)/2 213]; % This is the singlet (3.20 ppm, tcho4 in the paper); multiple is tcho3: [1310 1310 (243+191)/2 178];
relax.PCr   = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting; same as Cr
relax.PE    = [1340 1190 (119+86)/2 158];
relax.Scy   = [1340 1190 (125+107)/2 170];
relax.Tau   = [1340 1190 (123+102)/2 (123+102)/2]; % No WM estimate available; take GM estimate

% Check if metabolite name is in the look-up table
if isfield(relax, metName)
    T1_GM = relax.(metName)(1) * 1e-3;
    T1_WM = relax.(metName)(2) * 1e-3;
    T2_GM = relax.(metName)(3) * 1e-3;
    T2_WM = relax.(metName)(4) * 1e-3;
else
    % If not, use an average
    T1_GM = 1340 * 1e-3;
    T1_WM = 1190 * 1e-3;
    T2_GM = 140 * 1e-3;
    T2_WM = 169 * 1e-3;
end

end

%%% / Lookup function for metabolite relaxation times %%%

%%% Function to create metabolite overview in MATLAB table format %%%
function [MRSCont] = osp_createTable(MRSCont, qtfyType, getResults)
    if ~(strcmp(qtfyType, 'AlphaCorrWaterScaled') || strcmp(qtfyType, 'AlphaCorrWaterScaledGroupNormed'))
        if ~strcmp(qtfyType, 'amplMets')
            % Extract metabolite names from basisset
            for ll = 1:length(getResults)
            names = MRSCont.quantify.metabs.(getResults{ll});

            conc = zeros(MRSCont.nDatasets,length(names));

                for kk = 1:MRSCont.nDatasets
                    conc(kk,:) = MRSCont.quantify.(getResults{ll}).(qtfyType){kk};
                end
                % Save back to Osprey data container
                if isfield(MRSCont, 'exclude')
                    if~isempty(MRSCont.exclude)
                        conc(MRSCont.exclude,:) = [];
                    end
                end
                MRSCont.quantify.tables.(getResults{ll}).(qtfyType)  = array2table(conc,'VariableNames',names);
            end
        else
            % Extract metabolite names from basisset
            for ll = 1:length(getResults)
            names = MRSCont.quantify.metabs.(getResults{ll});

            conc = zeros(MRSCont.nDatasets,length(names));

                for kk = 1:MRSCont.nDatasets
                    conc(kk,:) = MRSCont.quantify.(qtfyType){kk}.(getResults{ll});
                end
                % Save back to Osprey data container
                if isfield(MRSCont, 'exclude')
                    if~isempty(MRSCont.exclude)
                        conc(MRSCont.exclude,:) = [];
                    end
                end
                MRSCont.quantify.tables.(getResults{ll}).(qtfyType)  = array2table(conc,'VariableNames',names);
            end            
        end
    else
            % Extract metabolite names from basisset
        names = {'GABA'};
        if ~strcmp(MRSCont.opts.fit.coMM3, 'none')
            names = {'GABA','GABAplus'};    
        end


        conc = zeros(MRSCont.nDatasets,length(names));
        for ll = 1:length(getResults)
            for kk = 1:MRSCont.nDatasets
                conc(kk,:) = MRSCont.quantify.(getResults{ll}).(qtfyType){kk};
            end
            % Save back to Osprey data container
            if isfield(MRSCont, 'exclude')
                if~isempty(MRSCont.exclude)
                    conc(MRSCont.exclude,:) = [];
                end
            end
            MRSCont.quantify.tables.(getResults{ll}).(qtfyType)  = array2table(conc,'VariableNames',names);
        end
    end
end
%%% Function to create metaboite overview in MATLAB table format %%%
