function [MRSCont] = OspreyQuantifyMRSI(MRSCont,MetabSpecName)
%% [MRSCont] = OspreyQuantifyMRSI(MRSCont)
%   This function transforms the raw amplitude parameters determined during
%   OspreyFit into MRSI  maps.
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
%       Dr. Helge Zollner (Johns Hopkins University, 2021-01-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%

outputFolder = MRSCont.outputFolder;
% diary(fullfile(outputFolder, 'LogFile.txt'));
% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    msg = 'Trying to quantify data, but data have not been modelled yet. Run OspreyFit first.';
    fprintf(msg);
    error(msg);    
end


%%% 0. CHECK WHICH QUANTIFICATIONS CAN BE DONE %%%
if isfield(MRSCont.opts, 'MRSI') && isfield(MRSCont.opts.MRSI,'WaterVisibility')
    WaterVisibility = MRSCont.opts.MRSI.WaterVisibility;
else
    WaterVisibility = 0.65;
end

% tCr ratios can always be calculated (unless the unlikely case that Cr is
% not in the basis set, a case we'll omit for now).
qtfyCr = 1;

% Check which types of metabolite data are available
if MRSCont.flags.isUnEdited
    getResults = {'metab'};
    tCrNorm = {'metab'};
elseif MRSCont.flags.isMEGA
    getResults = {'metab', 'metab'};
    tCrNorm = {'metab'};
elseif MRSCont.flags.isHERMES
    getResults = {'metab', 'metab', 'metab'};
    tCrNorm = {'metab'};
elseif MRSCont.flags.isHERCULES
    getResults = {'metab', 'metab', 'metab'};
    tCrNorm = {'metab'};
end

% Check which types of water data are available
if [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [1 1]
        % If both water reference and short-TE water data have been
        % provided, use the one with shorter echo time.
        qtfyH2O     = 1;
        getResultsWater = {'water'};
        waterType = 'w';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [1 0]
        % If only one type of water data has been provided, use it.
        qtfyH2O     = 1;
        getResultsWater = {'ref'};
        waterType = 'ref';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [0 1]
        % If only one type of water data has been provided, use it.
        qtfyH2O     = 1;
        getResultsWater = {'water'};
        waterType = 'w';
elseif [MRSCont.flags.hasRef MRSCont.flags.hasWater] == [0 0]
        % If no water ref has been provided, only tCr ratios can be
        % provided.
        qtfyH2O     = 0;
end

% Get the fieldstrength for proper relaxation correction
if qtfyH2O
    Bo = MRSCont.raw{1}.Bo;  
    if (Bo >= 2.8 && Bo < 3.1)
            Bo = '3T';
    else
            Bo = '7T';
    end
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
saveDestination = fullfile(outputFolder,'nii-export',['fit_raw_' MetabSpecName]);
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end




%% Loop over all datasets
QuantifyTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyQuant',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
    
    % Get all the model results from fit object
    ModelMatrix = MRSCont.fit.(getResults{1});   
    ModelMatrix=flip(ModelMatrix,1);
    if (MRSCont.raw{1}.nZvoxels > 1)
        ModelMatrix = flip(ModelMatrix,length(size(ModelMatrix)));
    end    
   MRSI_model = get_MRSI_amplitudes(ModelMatrix); 
    non_zero = find(~cellfun('isempty', ModelMatrix));
    temp = ModelMatrix{non_zero(1)};
    % metsName = temp.BasisSets.names(logical(temp.BasisSets.includeInFit(end,:)));
    metsName = temp.Model{end}.CRLB.Properties.VariableNames;
    
    
    if qtfyH2O
        WaterModelMatrix = MRSCont.fit.(getResultsWater{1});
        WaterModelMatrix = flip(WaterModelMatrix ,1);
        if (MRSCont.raw_w{1}.nZvoxels > 1)
            WaterModelMatrix = flip(WaterModelMatrix,length(size(WaterModelMatrix)));
        end
        MRSI_model_water = get_MRSI_amplitudes(WaterModelMatrix);         
        MRSI_model_water.amplitudes = squeeze(sum(MRSI_model_water.amplitudes,1));             
    end
    if qtfyCr
        ModelMatrix2 = MRSCont.fit.(tCrNorm{1});        
        ModelMatrix2=flip(ModelMatrix2,2);
        if (MRSCont.raw{1}.nZvoxels > 1)
            ModelMatrix2 = flip(ModelMatrix2,length(size(ModelMatrix2)));
        end
        non_zero = find(~cellfun('isempty', ModelMatrix2));
        temp = ModelMatrix2{non_zero(1)};
        MRSI_model2 = get_MRSI_amplitudes(ModelMatrix2); 
        metab_namesNorm = temp.Model{1}.CRLB.Properties.VariableNames;
        tCrCombinationNames = {'tNAA','tCr','tCr_methyl_only'};
        for ll = 1 : 3
            Idx_1 = find(strcmp(metab_namesNorm,tCrCombinationNames{ll})); 
            if ~isempty(Idx_1)
                NormIndex = Idx_1;
            end
        end
    end
    
    if qtfyH2O
         % Get repetition times
        metsTR  = MRSCont.processed.A{kk}.tr * 1e-3;
        waterTR = MRSCont.processed.(waterType){kk}.tr * 1e-3;
        % Get echo times
        metsTE  = MRSCont.processed.A{kk}.te * 1e-3;
        waterTE = MRSCont.processed.(waterType){kk}.te * 1e-3;
         % Calculate factor for water-scaled, but not tissue-corrected metabolite levels
        rawWaterScaledFactor = quantH2O(metsName, metsTR, waterTR, metsTE, waterTE,Bo, WaterVisibility);
    end

    if qtfyCSF
        fCSF = squeeze(MRSCont.seg.tissue.fCSF(kk,:,:,:));
    end

    if qtfyTiss
        % Get repetition times
        metsTR  = MRSCont.processed.A{kk}.tr * 1e-3;
        waterTR = MRSCont.processed.(waterType){kk}.tr * 1e-3;
        % Get echo times
        metsTE  = MRSCont.processed.A{kk}.te * 1e-3;
        waterTE = MRSCont.processed.(waterType){kk}.te * 1e-3;
         % Calculate factor for water-scaled, but not tissue-corrected metabolite levels
        TissCorrWaterScaledFactor = quantTiss(metsName, metsTR, waterTR, metsTE, waterTE, squeeze(MRSCont.seg.tissue.fGM(kk,:,:,:)), squeeze(MRSCont.seg.tissue.fWM(kk,:,:,:)), squeeze(MRSCont.seg.tissue.fCSF(kk,:,:,:)),Bo);
    end

    if (MRSCont.raw{kk}.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
        reorder = flip(1:data.nZvoxels);
        for ll = 1 : data.nZvoxels
                ToExport = data;
                if isfield(ToExport.geometry,'slice_distance')
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                end
                mkdir(fullfile(saveDestination,['slice_' num2str(reorder(ll))],'fit'))
                out.hdr = ToExport.nii_mrs.hdr;
                out.hdr.dim(1) = 3;
                out.hdr.dim(2) = 1;
                out.hdr.pixdim(5) = 1;
                out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
                if ~isempty(WaterModelMatrix)
                    mkdir(fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','rawWaterScaled'))
                end
                mkdir(fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','amplitudes'))
                mkdir(fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','tCr'))
                mkdir(fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','CRLBs'))
                for mm = 1 : length(metsName)
                    if qtfyH2O
                        out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./MRSI_model_water.amplitudes(:,:,ll) *rawWaterScaledFactor(mm));         
                    
                        out.img(isnan(out.img)) =0;
                        out.img(isinf(out.img)) =0;
                        nii_tool('save', out, fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','rawWaterScaled',[metsName{mm}  '.nii.gz']));
                    end

                    out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./(squeeze(MRSI_model2.amplitudes(NormIndex,:,:,ll))));
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','tCr',[metsName{mm}  '.nii.gz']));
    
                    out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll)));
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','amplitudes',[metsName{mm}  '.nii.gz']));
    
    
                    out.img = squeeze(MRSI_model.relCRLBs(mm,:,:,ll));
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    nii_tool('save', out, fullfile(saveDestination,['slice_' num2str(reorder(ll))],'concs','CRLBs',[metsName{mm}  '_CRLBs.nii.gz']));
                end
                
    
            shift = shift - 1;
        end
    else
        ToExport = MRSCont.processed.A{kk};
        mkdir(fullfile(saveDestination,'fit'))
        out.hdr = ToExport.nii_mrs.hdr;
        out.hdr.dim(1) = 3;
        out.hdr.dim(2) = 1;
        out.hdr.pixdim(5) = 1;
        out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
        if qtfyH2O
            mkdir(fullfile(saveDestination,'concs','rawWaterScaled'))
            mkdir(fullfile(saveDestination,'concs','rawWaterScaled_QCfilt'))
            mkdir(fullfile(saveDestination,'concs','rawWaterScaled_QC'))
        end
        if qtfyCSF
            mkdir(fullfile(saveDestination,'concs','CSFWaterScaled'))
            mkdir(fullfile(saveDestination,'concs','CSFWaterScaled_QCfilt'))
            mkdir(fullfile(saveDestination,'concs','CSFWaterScaled_QC'))
        end
        if qtfyTiss
            mkdir(fullfile(saveDestination,'concs','TissCorrWaterScaled'))
            mkdir(fullfile(saveDestination,'concs','TissCorrWaterScaled_QCfilt'))
            mkdir(fullfile(saveDestination,'concs','TissCorrWaterScaled_QC'))
        end
        mkdir(fullfile(saveDestination,'concs','QC'))
        mkdir(fullfile(saveDestination,'concs','amplitudes'))
        mkdir(fullfile(saveDestination,'concs','tCr'))
        mkdir(fullfile(saveDestination,'concs','amplitudes_QCfilt'))
        mkdir(fullfile(saveDestination,'concs','tCr_QCfilt'))
        mkdir(fullfile(saveDestination,'concs','amplitudes_QC'))
        mkdir(fullfile(saveDestination,'concs','tCr_QC'))
        mkdir(fullfile(saveDestination,'concs','CRLBs')) 

         % Export global QC maps
        if MRSCont.flags.didSeg
            brain_mask = squeeze(MRSCont.seg.tissue.brain(1,:,:,:));
            % brain_mask = flip(brain_mask,3);
        else
            if isfield(MRSCont.opts.MRSI,'MRSImask')
                brain_mask = MRSCont.opts.MRSI.MRSImask;
            else
                brain_mask = ones(MRSCont.raw{1}.nXvoxels,MRSCont.raw{1}.nYvoxels,MRSCont.raw{1}.nZvoxels);
            end
        end
        
        GlobalQC = brain_mask * 3;
        FWHM = MRSCont.quickMaps.A.FWHM;
        SNR = MRSCont.quickMaps.A.SNR;
        FWHM = flip(FWHM,1);
        SNR = flip(SNR,1);
        FWHM = flip(FWHM,3);
        SNR = flip(SNR,3);
        GlobalQC(FWHM>MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold)=1;
        GlobalQC((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & (SNR<MRSCont.opts.MRSI.Quantify.QC.SNRThreshold))=2;
        out.img = GlobalQC .* brain_mask;
        out.img(isnan(out.img)) =0;
        out.img(isinf(out.img)) =0;
        nii_tool('save', out, fullfile(saveDestination,'concs','QC',[  'GlobalQC_FWHM_SNR.nii.gz']));

        for mm = 1 : length(metsName)

            % Prepare metabolite QC maps
            CRLB = squeeze(MRSI_model.relCRLBs(mm,:,:,:));
            MetaboliteQC = brain_mask * 5;
            MetaboliteQC(FWHM>MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold)=1;
            MetaboliteQC((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                         (SNR<MRSCont.opts.MRSI.Quantify.QC.SNRThreshold))=2;
            MetaboliteQC((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                         (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                         (CRLB>MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold))=3;
            MetaboliteQC_FWHM_SNR_CRLB = MetaboliteQC;

            % Export raw amplitudes & QC map
            out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:)));        
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            percentile_val = prctile(out.img(:), MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold);
            nii_tool('save', out, fullfile(saveDestination,'concs','amplitudes',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.amplitudes.(metsName{mm}) = out.img;
            
            MetaboliteQC_FWHM_SNR_CRLB((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                                       (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                                       (CRLB<MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold) & ...
                                       (out.img>percentile_val))=4;

            out.img(out.img>percentile_val) =0;
            out.img = out.img  .* brain_mask;
            nii_tool('save', out, fullfile(saveDestination,'concs','amplitudes_QCfilt',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.amplitudes_QCfilt.(metsName{mm}) = out.img;

            out.img = MetaboliteQC_FWHM_SNR_CRLB;
            nii_tool('save', out, fullfile(saveDestination,'concs','amplitudes_QC',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.amplitudes_QC.(metsName{mm}) = out.img;
            

            % Export tCr ratios & QC map
            MetaboliteQC_FWHM_SNR_CRLB = MetaboliteQC;
            out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./(squeeze(MRSI_model2.amplitudes(end-1,:,:,:))));
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            percentile_val = prctile(out.img(:), MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold);
            nii_tool('save', out, fullfile(saveDestination,'concs','tCr',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.tCr.(metsName{mm}) = out.img;

            MetaboliteQC_FWHM_SNR_CRLB((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                                       (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                                       (CRLB<MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold) & ...
                                       (out.img>percentile_val))=4;

            out.img(out.img>percentile_val) =0;
            nii_tool('save', out, fullfile(saveDestination,'concs','tCr_QCfilt',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.tCr_QCfilt.(metsName{mm}) = out.img;

            out.img = MetaboliteQC_FWHM_SNR_CRLB;
            out.img = out.img  .* brain_mask;
            nii_tool('save', out, fullfile(saveDestination,'concs','tCr_QC',[metsName{mm}  '.nii.gz']));
            MRSCont.quantify.tCr_QC.(metsName{mm}) = out.img;
            
            % Export raw water scaled & QC map
            if qtfyH2O        
                MetaboliteQC_FWHM_SNR_CRLB = MetaboliteQC;
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./MRSI_model_water.amplitudes(:,:,:)) * rawWaterScaledFactor.(metsName{mm})(1);
                out.img(isnan(out.img)) =0;
                out.img(isinf(out.img)) =0;
                percentile_val = prctile(out.img(:), MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold);
                nii_tool('save', out, fullfile(saveDestination,'concs','rawWaterScaled',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.rawWaterScaled.(metsName{mm}) = out.img;

                MetaboliteQC_FWHM_SNR_CRLB((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                                           (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                                           (CRLB<MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold) & ...
                                           (out.img>percentile_val))=4;
    
                out.img(out.img>percentile_val) =0;
                nii_tool('save', out, fullfile(saveDestination,'concs','rawWaterScaled_QCfilt',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.rawWaterScaled_QCfilt.(metsName{mm}) = out.img;
    
                out.img = MetaboliteQC_FWHM_SNR_CRLB;
                out.img = out.img  .* brain_mask;
                nii_tool('save', out, fullfile(saveDestination,'concs','rawWaterScaled_QC',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.rawWaterScaled_QC.(metsName{mm}) = out.img;
            end
            
            % Export csf-corrected raw water scaled & QC map
            if qtfyCSF      
                MetaboliteQC_FWHM_SNR_CRLB = MetaboliteQC;
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./MRSI_model_water.amplitudes(:,:,:)) * rawWaterScaledFactor.(metsName{mm})(1) ./ (1-fCSF);
                out.img(isnan(out.img)) =0;
                out.img(isinf(out.img)) =0;
                percentile_val = prctile(out.img(:), MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold);
                nii_tool('save', out, fullfile(saveDestination,'concs','CSFWaterScaled',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.CSFrawWaterScaled.(metsName{mm}) = out.img;

                MetaboliteQC_FWHM_SNR_CRLB((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                                           (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                                           (CRLB<MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold) & ...
                                           (out.img>percentile_val))=4;
    
                out.img(out.img>percentile_val) =0;
                nii_tool('save', out, fullfile(saveDestination,'concs','CSFWaterScaled_QCfilt',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.CSFrawWaterScaled_QCfilt.(metsName{mm}) = out.img;
    
                out.img = MetaboliteQC_FWHM_SNR_CRLB;
                out.img = out.img  .* brain_mask;
                nii_tool('save', out, fullfile(saveDestination,'concs','CSFWaterScaled_QC',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.CSFrawWaterScaled_QC.(metsName{mm}) = out.img;
                
            end

            % Export tissue-corrected raw water scaled & QC map
            if qtfyTiss  
                MetaboliteQC_FWHM_SNR_CRLB = MetaboliteQC;
                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,:))./MRSI_model_water.amplitudes(:,:,:)) .* TissCorrWaterScaledFactor.(metsName{mm})(:,:,:);               
                out.img(isnan(out.img)) =0;
                out.img(isinf(out.img)) =0;     
                % out.img = nonlocalMeansDenoise(out.img);
                % out.img(isnan(out.img)) =0;
                % out.img(isinf(out.img)) =0; 
                percentile_val = prctile(out.img(:), MRSCont.opts.MRSI.Quantify.QC.PercentileThreshold);
                nii_tool('save', out, fullfile(saveDestination,'concs','TissCorrWaterScaled',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.TissCorrWaterScaled.(metsName{mm}) = out.img;

                MetaboliteQC_FWHM_SNR_CRLB((FWHM<MRSCont.opts.MRSI.Quantify.QC.FWHMThreshold) & ...
                                           (SNR>MRSCont.opts.MRSI.Quantify.QC.SNRThreshold) & ...
                                           (CRLB<MRSCont.opts.MRSI.Quantify.QC.CRLBThreshold) & ...
                                           (out.img>percentile_val))=4;
    
                out.img(out.img>percentile_val) =0;
                nii_tool('save', out, fullfile(saveDestination,'concs','TissCorrWaterScaled_QCfilt',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.TissCorrWaterScaled_QCfilt.(metsName{mm}) = out.img;
    
                out.img = MetaboliteQC_FWHM_SNR_CRLB;
                out.img = out.img  .* brain_mask;
                nii_tool('save', out, fullfile(saveDestination,'concs','TissCorrWaterScaled_QC',[metsName{mm}  '.nii.gz']));
                MRSCont.quantify.TissCorrWaterScaled_QC.(metsName{mm}) = out.img;
            end

            % Export CRLBs
            out.img = squeeze(MRSI_model.relCRLBs(mm,:,:,:));
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            nii_tool('save', out, fullfile(saveDestination,'concs','CRLBs',[metsName{mm}  '_CRLBs.nii.gz']));
            MRSCont.quantify.CRLBs.(metsName{mm}) = out.img;
           

        end
        if qtfyH2O
            out.img = squeeze(squeeze(MRSI_model_water.amplitudes(:,:,:)));
            out.img(isnan(out.img)) =0;
            out.img(isinf(out.img)) =0;
            nii_tool('save', out, fullfile(saveDestination,'concs','amplitudes',[  'water.nii.gz']));
        end
         MRSCont.quantify.GlobalQC =  GlobalQC .* brain_mask;
    end

    if qtfyH2O
        temp_img = squeeze(squeeze(MRSI_model_water.amplitudes(:,:,:)));
        temp_img(isnan(temp_img)) =0;
        temp_img(isinf(temp_img)) =0;
        MRSCont.quantify.water = temp_img;
    end

    % copy the correct fileTree logic for FSLeyes
    if qtfyTiss 
        copyfile(which(fullfile('mrsi','osprey_mrsi_qtfyTiss.tree')), saveDestination)
    elseif qtfyH2O
        copyfile(which(fullfile('mrsi','osprey_mrsi_qtfyH2O.tree')), saveDestination)
    else
        copyfile(which(fullfile('mrsi','osprey_mrsi_qtfyCr.tree')), saveDestination)
    end
    % Copy colorscheme
    copyfile(which(fullfile('mrsi','osprey_colourscheme.json')), saveDestination)
    

end
time = toc(QuantifyTime);
if MRSCont.flags.isGUI && isfield(progressText,'String')      
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf('... done.\n Elapsed time %f seconds\n',time);
MRSCont.runtime.Quantify = time;
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


%%% Calculate correction for water-scaled estimates %%%
function rawWaterScaledFactor = quantH2O(metsName, metsTR, waterTR, metsTE, waterTE,Bo,WaterVisibility)

% Define constants
PureWaterConc       = 1;            % mmol/L
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
for kk = 1:length(metsName)
    [T1_Metab_GM(kk), T1_Metab_WM(kk), T2_Metab_GM(kk), T2_Metab_WM(kk)] = lookUpRelaxTimes(metsName{kk},Bo);
    % average across GM and WM
    T1_Metab(kk) = mean([T1_Metab_GM(kk) T1_Metab_WM(kk)]);
    T2_Metab(kk) = mean([T2_Metab_GM(kk) T2_Metab_WM(kk)]);
    T1_Factor(kk) = (1-exp(-waterTR./T1_Water)) ./ (1-exp(-metsTR./T1_Metab(kk)));
    T2_Factor(kk) = exp(-waterTE./T2_Water) ./ exp(-metsTE./T2_Metab(kk));

    % Calculate
    rawWaterScaledFactor.(metsName{kk}) = PureWaterConc ...
        .* WaterVisibility .* T1_Factor(kk) .* T2_Factor(kk);
end

end
%%% /Calculate raw water-scaled estimates %%%



%%% Calculate tissue-corrected water-scaled estimates %%%
function TissCorrWaterScaledFactor = quantTiss(metsName,metsTR, waterTR, metsTE, waterTE, fGM, fWM, fCSF,Bo)
% This function calculates water-scaled, tissue-corrected metabolite
% estimates in molal units, according to Gasparovic et al, Magn Reson Med
% 55:1219-26 (2006).

% Define Constants
switch Bo
    case '3T'
        % Water relaxation 3 T
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
    case '7T'
        % Water relaxation 7 T
         % T1 from Rooney et al. 2007 (MRM)
         % T2 from Bartha et al. 2002 (MRM)
        T1w_WM    = 1.220;
        T2w_WM    = 0.050;
        T1w_GM    = 2.130;
        T2w_GM    = 0.055;
        T1w_CSF   = 4.425;
        T2w_CSF   = 0.141;
end

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

metsTE              = metsTE * 1e-3;    % convert to s
waterTE             = waterTE * 1e-3;   % convert to s
metsTR              = metsTR * 1e-3;    % convert to s
waterTR             = waterTR * 1e-3;   % convert to s


concW_GM    = 43.30*1e3;
concW_WM    = 36.08*1e3;
concW_CSF   = 53.84*1e3;
molal_concW = 55.51*1e3;

% Gasparovic et al. method
% Calculate molal fractions from volume fractions (equivalent to eqs. 5-7 in Gasparovic et al., 2006)
molal_fGM  = (fGM*concW_GM) ./ (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
molal_fWM  = (fWM*concW_WM) ./ (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
molal_fCSF = (fCSF*concW_CSF) ./ (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);

% Metabolites
for kk = 1:length(metsName)
    [T1_Metab_GM(kk), T1_Metab_WM(kk), T2_Metab_GM(kk), T2_Metab_WM(kk)] = lookUpRelaxTimes(metsName{kk},Bo);
    % average across GM and WM
    T1_Metab(kk) = mean([T1_Metab_GM(kk) T1_Metab_WM(kk)]);
    T2_Metab(kk) = mean([T2_Metab_GM(kk) T2_Metab_WM(kk)]);

    % Calculate water-scaled, tissue-corrected molal concentration
    % estimates
    TissCorrWaterScaledFactor.(metsName{kk})  = molal_concW ...
        .* (molal_fGM  * (1 - exp(-waterTR/T1w_GM)) * exp(-waterTE/T2w_GM) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk))) + ...
            molal_fWM  * (1 - exp(-waterTR/T1w_WM)) * exp(-waterTE/T2w_WM) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk))) + ...
            molal_fCSF * (1 - exp(-waterTR/T1w_CSF)) * exp(-waterTE/T2w_CSF) / ((1 - exp(-metsTR/T1_Metab(kk))) * exp(-metsTE/T2_Metab(kk)))) ./ ...
            (1 - molal_fCSF);
end

end
%%% /Calculate CSF-corrected water-scaled estimates %%%



%%% Lookup function for metabolite relaxation times %%%
function [T1_GM, T1_WM, T2_GM, T2_WM] = lookUpRelaxTimes(metName,Bo)

% Look up table below
switch Bo
    case '3T'
        % T1 values for NAA, Glu, Cr, Cho, Ins from Mlynarik et al, NMR Biomed
        % 14:325-331 (2001)
        % T1 for GABA from Puts et al, J Magn Reson Imaging 37:999-1003 (2013)
        % T2 values from Wyss et al, Magn Reson Med 80:452-461 (2018)
        % T2 values are averaged between OCC and pACC for GM; and PVWM for WM
        % Structure is [T1_GM T1_WM T2_GM T2_WM]
        relax.Asc   = [1340 1190 (125+105)/2 172];
        relax.Asp   = [1340 1190 (111+90)/2 148];
        relax.Cr    = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting
        relax.Cr_methyl_only    = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 
        relax.GABA  = [1310 1310 (102+75)/2 (102+75)/2]; % No WM estimate available; take GM estimate; both in good accordance with 88 ms reported by Edden et al
        relax.Glc   = [1340 1190 (117+88)/2 155]; % Glc1: [1310 1310 (128+90)/2 156];
        relax.Gln   = [1340 1190 (122+99)/2 168];
        relax.Glu   = [1270 1170 (135+122)/2 124];
        relax.Gly   = [1340 1190 (102+81)/2 152];
        relax.GPC   = [1300 1080 (274+222)/2 218]; % This is the Choline singlet (3.21 ppm, tcho2 in the paper); glycerol is tcho: [1310 1310 (257+213)/2 182]; % choline multiplet is tcho1: [1310 1310 (242+190)/2 178];
        relax.GPC_pCh2_only   = [1300 1080 (274+222)/2 218]; % This is the Choline singlet (3.21 ppm, tcho2 in the paper); 
        relax.GSH   = [1340 1190 (100+77)/2 145]; % This is the cysteine signal (GSH1 in the paper), glycine is GSH: [1310 1310 (99+72)/2 145]; % glutamate is GSH2: [1310 1310 (102+76)/2 165];
        relax.Lac   = [1340 1190 (110+99)/2 159];
        relax.Ins   = [1230 1010 (244+229)/2 161];
        relax.NAA   = [1470 1350 (253+263)/2 343]; % This is the 2.008 ppm acetyl signal (naa in the paper); aspartyl is naa1: [1310 1310 (223+229)/2 310];
        relax.NAA_Acetyl_only   = [1470 1350 (253+263)/2 343]; % This is the 2.008 ppm acetyl signal (naa in the paper); 
        relax.NAAG  = [1340 1190 (128+107)/2 185]; % This is the 2.042 ppm acetyl signal (naag in the paper); aspartyl is naag1: [1310 1310 (108+87)/2 180]; % glutamate is NAAG2: [1310 1310 (110+78)/2 157];
        relax.NAAG_Acetyl_only  = [1340 1190 (128+107)/2 185]; % This is the 2.042 ppm acetyl signal (naag in the paper); aspartyl is naag1: [1310 1310 (108+87)/2 180]; 
        relax.PCh   = [1300 1080 (274+221)/2 213]; % This is the singlet (3.20 ppm, tcho4 in the paper); multiple is tcho3: [1310 1310 (243+191)/2 178];
        relax.PCh_trimethyl_only   = [1300 1080 (274+221)/2 213]; % This is the singlet (3.20 ppm, tcho4 in the paper); 
        relax.PCr   = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting; same as Cr
        relax.PCr_ch3_only   = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 
        relax.PE    = [1340 1190 (119+86)/2 158];
        relax.Scy   = [1340 1190 (125+107)/2 170];
        relax.Tau   = [1340 1190 (123+102)/2 (123+102)/2]; % No WM estimate available; take GM estimate
        relax.tNAA  = [(1470+1340)/2 (1350+1190)/2 (253+263+128+107)/4 (343+185)/2]; % Mean values from NAA + NAAG
        relax.tNAA_Acetyl_only  = [(1470+1340)/2 (1350+1190)/2 (253+263+128+107)/4 (343+185)/2]; % Mean values from NAA + NAAG
        relax.tCr  = [(1460+1460)/2 (1240+1240)/2 (148+144+148+144)/4 (166+166)/2]; % Mean values from Cr + PCr
        relax.tCr_methyl_only  = [(1460+1460)/2 (1240+1240)/2 (148+144+148+144)/4 (166+166)/2]; % Mean values from Cr + PCr
        relax.tCho  = [(1300+1080)/2 (1080+1080)/2 (274+222+274+221)/4 (218+213)/2]; % Mean values from GPC + PCh
        relax.tCho_pCh2_only  = [(1300+1080)/2 (1080+1080)/2 (274+222+274+221)/4 (218+213)/2]; % Mean values from GPC + PCh
        relax.Glx  = [(1340+1270)/2 (1190+1170)/2 (122+99+135+122)/4 (168+124)/2]; % Mean values from Glu + Glx

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
    case '7T'
         % T2 values of water, NAA, tCr, tCho, Scyllo, Ins, Glu,GSH, Ins,
         % and Tau are taken from Marjanska et al. 2011 (NMR
         % 10.1002/nbm.1754). It was averaged across 4 regions OCC, SM1, BG, CER 
         % Penner et al (2014) https://doi.org/10.1002/mrm.25380 for 
        relax.Asc   = [1530 1484 127 128]; % This is the average from tNAA, tCr, tCho, Glx, and Ins
        relax.Asp   = [1530 1484 127 128]; % This is the average from tNAA, tCr, tCho, Glx, and Ins
        relax.Cr    = [1740 1780 107 107]; % Taken from tCr
        relax.GABA  = [1334 1334 87 87]; % Andreychenko et al. (2012) 10.1002/nbm.2997
        relax.Glc   = [1530 1484 127 128]; % This is the average from tNAA, tCr, tCho, Glx, and Ins 
        relax.Gln   = [1640 1740 107 107]; %T1 from Mlynarik et al. (2012) 10.1002/mrm.24352 % T2 as Gln
        relax.Glu   = [1610 1750 107 117]; % T1 from Mlynarik et al. (2012) 10.1002/mrm.24352, T2 from https://doi.org/10.1371/journal.pone.0215210
        relax.Gly   = [1530 1484 127 128]; % This is the average from tNAA, tCr, tCho, Glx, and Ins
        relax.GPC   = [1510 1320  153 153]; % Taken from tCho
        relax.GSH   = [1140 1060 79 79]; % Entire molecule; T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
        relax.Lac   = [1530 1484 182 182]; % The use of MEGA-sLASER with J-refocusing echo time extension to measure the proton T2 of lactate in healthy human brain at 7 T ISMRM
        relax.Ins   = [1280 1190 111 111]; %T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
        relax.NAA   = [1780 1830 155 155]; % This is the 2.008 ppm acetyl signal (naa in the paper); aspartyl is naa1: [1310 1310 110 110]; T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
        relax.NAAG  = [1210 940 155 155]; % This is the 2.042 ppm acetyl signal (naag in the paper); aspartyl is naag1: [1310 1310 (108+87)/2 180]; % glutamate is NAAG2: [1310 1310 (110+78)/2 157]; 
        relax.PCh   = [1510 1320  153 153]; % Taken from tCho
        relax.PCr   = [1740 1780 107 107]; % Taken from tCr
        relax.PE    = [1310 1320]; %T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
        relax.Scy   = [1310 1230 105 105]; %T1 from Mlynarik et al. (2012) 10.1002/mrm.24352 T2 from https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.24352
        relax.Tau   = [2150 2090 97 97]; %T1 from Mlynarik et al. (2012) 10.1002/mrm.24352 % T2 as Gln
        relax.tNAA  = [1495 1385 155 155]; % Mean values from NAA + NAAG
        relax.tCr  = [1740 1780 107 107]; % The singlet peak ar 3 ppm. 3.9 ppm peak values are [1240 1190 94 94] %T1 from Mlynarik et al. (2012)
        relax.tCho  = [1510 1320  153 153]; % Entire molecule; T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
        relax.Glx  = [1625 1745 107 112]; % Mean values from Glu + Glx

        % Check if metabolite name is in the look-up table
        if isfield(relax, metName)
            T1_GM = relax.(metName)(1) * 1e-3;
            T1_WM = relax.(metName)(2) * 1e-3;
            T2_GM = relax.(metName)(3) * 1e-3;
            T2_WM = relax.(metName)(4) * 1e-3;
        else
            % If not, use an average
            T1_GM = 1530 * 1e-3; % This is the average from tNAA, tCr, tCho, Glx, and Ins
            T1_WM = 1484 * 1e-3; % This is the average from tNAA, tCr, tCho, Glx, and Ins
            T2_GM = 127 * 1e-3; % This is the average from tNAA, tCr, tCho, Glx, and Ins
            T2_WM = 128 * 1e-3; % This is the average from tNAA, tCr, tCho, Glx, and Ins
        end
end

end

%%% / Lookup function for metabolite relaxation times %%%

