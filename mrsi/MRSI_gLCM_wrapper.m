function [MRSCont] = MRSI_gLCM_wrapper(MRSCont, ModelProcedureFileMetabolites, ModelProcedureFileWater, BasisSetFile, MetabSpecName, parallel_flag,zero_fill)
%% Wrapper for parallel processing of MRSI data in the MRSContainer
if isfield(MRSCont,'MRSImask')
    mask = MRSCont.MRSImask;
else
    mask = ones(MRSCont.processed.(MetabSpecName){1}.nXvoxels,MRSCont.processed.(MetabSpecName){1}.nYvoxels,MRSCont.processed.(MetabSpecName){1}.nZvoxels);
end
if zero_fill
    no_zero_fill = 0;
else
    no_zero_fill = 1;
end
MRSCont.runtime.Fit = 0;
outputFolder = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if parallel_flag
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj); % Shut down and delete the parallel pool if it exists
        disp('Old parallel pool shut down successfully.');
    end
    p = parpool(); % Start parpool once per session
end
%% Water model
if MRSCont.flags.hasWater
    waterFitTime = tic;
    MRSCont.processed.w{1}.nucleus = {'1H'};
    x_vec=[];
    y_vec=[];
    z_vec=[];

    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.w{1}.nZvoxels
        for x = 1 : MRSCont.processed.w{1}.nXvoxels
            for y = 1 : MRSCont.processed.w{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.w{1},[x y z]);               
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    WaitMessage = parfor_wait(length(data_vec),'Waitbar',true,'ReportInterval',round(length(data_vec)/20),'CurrentStep','Get water FWHM inital');
    if parallel_flag
        % p = parpool();
        parfor vx = 1:  length(data_vec)
            if mask_vec(vx) >= 1
                [~, refFWHM] = osp_XReferencing(data_vec{vx},4.68,1,[0 9.36],0);
                data_vec{vx}.FWHM = refFWHM * data_vec{vx}.txfrq*1e-6;
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                                'Manufacturer','OriginalFile','software','nii_mrs'};
                    data_vec{vx} = rmfield(data_vec{vx},fields);
            end
            WaitMessage.Send;
        end
        % delete(p);
    else
        for vx = 1:  length(data_vec)
            if mask_vec(vx) >= 1
                [~, refFWHM] = osp_XReferencing(data_vec{vx},4.68,1,[0 9.36],0);
                data_vec{vx}.FWHM = refFWHM * data_vec{vx}.txfrq*1e-6;
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                                'Manufacturer','OriginalFile','software','nii_mrs'};
                    data_vec{vx} = rmfield(data_vec{vx},fields);
            end
            WaitMessage.Send;
        end
    end
    WaitMessage.Destroy;


    % Create basis with matching dwelltime 
    load(BasisSetFile{1});                                        % Assume it is the first one ...

    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    if zero_fill
        DataToModel = op_zeropad(data_vec{1},2);   
    else
        DataToModel = data_vec{1};
    end
    BASIS = fit_resampleBasis(DataToModel, BASIS); 
    BASIS.centerFreq = 4.68;


    % Get data scale
    scaleData = max(real(MRSCont.processed.(MetabSpecName){1}.specs(MRSCont.processed.(MetabSpecName){1}.ppm > -2 & MRSCont.processed.(MetabSpecName){1}.ppm < 10 ,:,:)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    % scaleData=0;

    try
        BASIS = rmfield(BASIS,'specs');
    catch
    end
    ToKeep = find(strcmp(BASIS.name,'H2O'));
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.fids(:,end+1)= BASIS.fids(:,1);
    BASIS.fids(:,end+1)= BASIS.fids(:,1);
    BASIS.fids(:,end+1)= BASIS.fids(:,1);
    BASIS.name{end+1}= 'H2O_A';
    BASIS.name{end+1}= 'H2O_B';
    BASIS.name{end+1}= 'H2O_C';
    BASIS.nMets = 4;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);


    % Model water data 
    ModelProced = ModelProcedureFileWater;

    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end


    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    mask_vec(mask_vec==0)=[];


    ModelProcedureCell = cell(1,length(data_vec));
    tempFWHM = 2;
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            if ~isnan(data_vec{vx}.FWHM)
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = data_vec{vx}.FWHM;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ex = data_vec{vx}.FWHM;
            tempFWHM = data_vec{vx}.FWHM;
            else
                ModelProcedureCell{vx} = ModelProcedure;
                data_vec{vx} = rmfield(data_vec{vx},'FWHM');
                ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = tempFWHM;
                ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ex = tempFWHM;
            end
        end
    end


    water = cell(1,length(data_vec));
    final_water = cell(1,length(data_vec));
    water_shifts = zeros(1,length(data_vec));
    WaitMessage = parfor_wait(length(data_vec),'Waitbar',true,'ReportInterval',round(length(data_vec)/20),'CurrentStep','Model Water');
    if parallel_flag
        % p = parpool();
        D = parallel.pool.Constant(data_vec);
        M = parallel.pool.Constant(ModelProcedureCell);
        S = parallel.pool.Constant(scaleData);


        tstart = tic;
        parfor vx = 1:  length(data_vec)
                water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,no_zero_fill,S.Value,0,0,BASIS,1);
                WaitMessage.Send;
        end
        modeltime_water = toc(tstart);
        WaitMessage.Destroy;
        % delete(p);
        water_temp = cell(MRSCont.processed.w{1}.nXvoxels,MRSCont.processed.w{1}.nYvoxels,MRSCont.processed.w{1}.nZvoxels);
        for vx = 1:  length(data_vec)
           water{vx}.economizeStorage(1,1);                         % Remove basis set and jacobians
           water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = water{vx};
           [~,index] = max(water{vx}.Model{1,1}.parsOut.metAmpl);
           water_shifts(vx) = water{vx}.Model{1,1}.parsOut.freqShift(index);
        end
        water = water_temp;
    else
        tstart = tic;
        for vx = 1:  length(data_vec)
                water(vx) = Osprey_gLCM(data_vec(vx),ModelProcedureCell{vx},0,0,scaleData,0,0,BASIS,1);
                WaitMessage.Send;
        end
        WaitMessage.Destroy;
        modeltime_water = toc(tstart);
        water_temp = cell(MRSCont.processed.w{1}.nXvoxels,MRSCont.processed.w{1}.nYvoxels,MRSCont.processed.w{1}.nZvoxels);
        for vx = 1:  length(data_vec)
           water{vx}.economizeStorage(1,1);                         % Remove basis set and jacobians
           water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = water{vx};
           [~,index] = max(water{vx}.Model{1,1}.parsOut.metAmpl);
           water_shifts(vx) = water{vx}.Model{1,1}.parsOut.freqShift(index);
        end
        water = water_temp;
    end

    time = toc(waterFitTime);
    MRSCont.runtime.FitWater = time;
    MRSCont.runtime.Fit = time;
    MRSCont.fit.water = water;
    MRSCont.fit.scale = scaleData;
    export_model_to_NII(MRSCont,MRSCont.processed.w{1},MRSCont.fit.water ,fullfile(outputFolder,'nii-export','fit_raw_w'));
end

%% Metabolite model

    metFitTime = tic;
    MRSCont.processed.(MetabSpecName){1}.nucleus = {'1H'};
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.(MetabSpecName){1}.nZvoxels
        for x = 1 : MRSCont.processed.(MetabSpecName){1}.nXvoxels
            for y = 1 : MRSCont.processed.(MetabSpecName){1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.(MetabSpecName){1},[x y z]);
                if isfield('watersupp',data_vec{ind})
                    fields = {'nXvoxels','nYvoxels','nZvoxels','seq',...
                                'nii_mrs','specReg','watersupp',...
                                'refFWHM','refShift'};
                else
                    fields = {'nXvoxels','nYvoxels','nZvoxels','seq',...
                                'nii_mrs','specReg',...
                                'refFWHM','refShift'};
                end
                data_vec{ind} = rmfield(data_vec{ind},fields); 
                data_vec{ind}.centerFreq = data_vec{ind}.centerFreq(1);
                data_vec{ind}.txfrq = data_vec{ind}.txfrq(1);
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end    

    if zero_fill
        DataToModel = op_zeropad(data_vec{1},2);   
    else
        DataToModel = data_vec{1};
    end

    ModelProcedure = jsonToStruct(ModelProcedureFileMetabolites);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end

    model = cell(1,length(data_vec));
    ModelProcedureCell = cell(1,length(data_vec));

    for vx = 1:  length(data_vec)
        ModelProcedureCell{vx} = ModelProcedure;
    end


    if length(BasisSetFile) == 1
        load(BasisSetFile{1});
        BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
        BASIS = fit_sortBasisSet(BASIS);  % Sort according to Osprey standard
        BASIS = fit_resampleBasis(DataToModel, BASIS); 
        if isfield(ModelProcedure.basisset,'opts')                          % Apply options to basis, e.g. pick subspectra
            if isfield(ModelProcedure.basisset.opts,'index')               % Index for subspectra
                BASIS.fids = BASIS.fids(:,:,ModelProcedure.basisset.opts.index);  % Pick fids according to index
                BASIS.specs = BASIS.specs(:,:,ModelProcedure.basisset.opts.index); % Pick specs according to index
                BASIS.sz = size(BASIS.fids);                              % Recalculate size entry
                try
                    basisSet.nExtra = basisSet.sz(4);                               % Update extra dimension
                catch
                end
            end
        end
    else
        for bb = 1 : length(BasisSetFile)
            if bb == 1
                load(BasisSetFile{bb});          % Load basis set
                BASIS = fit_sortBasisSet(BASIS);  % Sort according to Osprey standard
                BASIS = fit_resampleBasis(DataToModel, BASIS); 
            else
                basisSetToAdd = load(BasisSetFile{bb});     % Load basis set
                basisSetToAdd = basisSetToAdd.BASIS;
                basisSetToAdd = fit_sortBasisSet(basisSetToAdd);  % Sort according to Osprey standard
                basisSetToAdd = fit_resampleBasis(DataToModel, basisSetToAdd); 
                BASIS.fids = cat(3,BASIS.fids,basisSetToAdd.fids);    % Concatenate time domain basis functions for 2D fit
                BASIS.specs = cat(3,BASIS.specs,basisSetToAdd.specs); % Concatenate frequency domain basis functions for 2D fit
            end
        end
    end
     BASIS.centerFreq = 4.68;


    if isfield(MRSCont,'fit') 
        if ~isfield(MRSCont.fit,'scale')
            scaleData = max(real(MRSCont.processed.(MetabSpecName){1}.specs(MRSCont.processed.(MetabSpecName){1}.ppm > -2 & MRSCont.processed.(MetabSpecName){1}.ppm < 10 ,:,:)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
        else
            scaleData =  MRSCont.fit.scale;
        end
    else
        scaleData = max(real(MRSCont.processed.(MetabSpecName){1}.specs(MRSCont.processed.(MetabSpecName){1}.ppm > -2 & MRSCont.processed.(MetabSpecName){1}.ppm < 10 ,:,:)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    end
    % scaleData=0;
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    mask_vec(mask_vec==0)=[];
    model = cell(1,length(data_vec));



    WaitMessage = parfor_wait(length(data_vec),'Waitbar',true,'ReportInterval',round(length(data_vec)/20),'CurrentStep','Model Metabolites');
    if parallel_flag
        % p = parpool();
        D = parallel.pool.Constant(data_vec);
        M = parallel.pool.Constant(ModelProcedureCell);
        S = parallel.pool.Constant(scaleData);
        tstart = tic;
        parfor vx = 1:  length(data_vec)
                model(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,no_zero_fill,S.Value,0,0,BASIS,1);
                WaitMessage.Send;
        end
        WaitMessage.Destroy;
        modeltime_final_none = toc(tstart);
        delete(p);
        model_temp = cell(MRSCont.processed.(MetabSpecName){1}.nXvoxels,MRSCont.processed.(MetabSpecName){1}.nYvoxels,MRSCont.processed.(MetabSpecName){1}.nZvoxels);
        for vx = 1:  length(data_vec)
           model{vx}.economizeStorage(1,1);                         % Remove basis set and jacobians
           model_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = model{vx};
        end
        model = model_temp;       
    else
        tstart = tic;
        for vx = 1:  length(data_vec)
                model(vx) = Osprey_gLCM(data_vec(vx),ModelProcedureCell{vx},0,no_zero_fill,scaleData,0,0,BASIS,1);
                WaitMessage.Send;
        end
        WaitMessage.Destroy;
        modeltime_final_none = toc(tstart);
        model_temp = cell(MRSCont.processed.(MetabSpecName){1}.nXvoxels,MRSCont.processed.(MetabSpecName){1}.nYvoxels,MRSCont.processed.(MetabSpecName){1}.nZvoxels);
        for vx = 1:  length(data_vec)
           model{vx}.economizeStorage(1,1);                         % Remove basis set and jacobians 
           model_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = model{vx};
        end
        model = model_temp;
    end

    MRSCont.fit.metab = model;
    MRSCont.flags.didFit = 1;
    time = toc(metFitTime);
    MRSCont.runtime.FitMet = time;

    MRSCont.runtime.Fit = MRSCont.runtime.Fit + MRSCont.runtime.FitMet;

    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    %% Export the model results
    export_model_to_NII(MRSCont,MRSCont.processed.(MetabSpecName){1},MRSCont.fit.metab ,fullfile(outputFolder,'nii-export',['fit_raw_metab']));
end
    