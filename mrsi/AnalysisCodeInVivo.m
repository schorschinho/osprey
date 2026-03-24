[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)
%%

for jobs = 19: length(matFiles)
    load(matFiles{jobs}); 
    load(maskFiles{jobs});

    MRSCont.processed.w{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    % mask = squeeze(mask(:,:,2));
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

    p = parpool('Dell');
    parfor vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            [~, refFWHM] = osp_XReferencing(data_vec{vx},4.68,1,[0 9.36],0);
            data_vec{vx}.FWHM = refFWHM * data_vec{vx}.txfrq*1e-6;
            fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs'};
                data_vec{vx} = rmfield(data_vec{vx},fields);
        end
    end
    delete(p);
    
    
    % Create basis with matching dwelltime 
    load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    DataToModel = op_zeropad(data_vec{1},2);   
    BASIS = fit_resampleBasis(DataToModel, BASIS); 
    BASIS.centerFreq = 4.68;
    
    
    % Get data scale
    % scaleData = max(real(MRSCont.processed.A{1}.specs(MRSCont.processed.A{1}.ppm > -2 & MRSCont.processed.A{1}.ppm < 10 ,:,:,2)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    % save(fullfile([folder], 'scaleData.mat'),'scaleData','-v7.3');
load(fullfile([folder], 'scaleData.mat'));

    BASIS = rmfield(BASIS,'specs');
    ToKeep = [16,29,30,31];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 4;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);

    clear MRSCont
    
    % Model water data 
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_water_in_vivo.json';
    
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
    

    ModelProcedureCell(mask_vec==0)=[];

    water = cell(1,length(data_vec));
    
    final_water = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_water = toc(tstart);
    delete(p);
    water_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = water{vx};
    end
    water = water_temp;
    save(fullfile([folder], 'water.mat'),'water','modeltime_water','-v7.3');
end
%%
% Model metabolite data with metab water prior knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)

for jobs = 4: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});                                  

    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'water.mat'))
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                water_vec(ind) = water(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end

    %Create basis with matching dwelltime 
    % load('/Users/helge/Documents/GitHub/osprey/fit/basissets/3T/philips/unedited/press/20/basis_philips_press20.mat');
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS.centerFreq = 4.68;
    BASIS = rmfield(BASIS,'specs');
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
    
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:18) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:18) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:18) = 1;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:18) = 0.5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:18) = 5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,19:26) = 6.44;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,19:26) = 6.44;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,19:26) = 2.12;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,19:26) = 0.5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,19:26) = 15;
    % 
    % ModelProcedure.Steps{1}.parametrizations.freqShift.init = 0;
    % ModelProcedure.Steps{1}.parametrizations.freqShift.lb = -10;
    % ModelProcedure.Steps{1}.parametrizations.freqShift.ub = 10;
    % ModelProcedure.Steps{1}.parametrizations.freqShift.ex = 0;
    % ModelProcedure.Steps{1}.parametrizations.freqShift.sd = 3;
    
    
    ModelProcedureCell = cell(1,length(data_vec));
    gauss_vec = ones(1,length(data_vec))*NaN;
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            gauss_vec(vx) = water_vec{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end
    

    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = water_vec{vx}.Model{1, 1}.parsOut.gaussLB;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = water_vec{vx}.Model{1, 1}.parsOut.gaussLB;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = water_vec{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    ModelProcedureCell(mask_vec==0)=[];
    mask_vec(mask_vec==0)=[];
    
    
    final_water = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_water = toc(tstart);
    final_water_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_water{vx};
    end
    final_water = final_water_temp;
    save(fullfile([folder], 'final_water.mat'),'final_water','modeltime_final_water','-v7.3');
    delete(p);
end
%%
% Model metabolite data with metab water prior knowledge and expectation
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});

    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'water.mat'))
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                water_vec(ind) = water(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end

    %Create basis with matching dwelltime 
    % load('/Users/helge/Documents/GitHub/osprey/fit/basissets/3T/philips/unedited/press/20/basis_philips_press20.mat');
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS.centerFreq = 4.68;
    BASIS = rmfield(BASIS,'specs');
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;

    
    ModelProcedureCell = cell(1,length(data_vec));
    gauss_vec = ones(1,length(data_vec))*NaN;
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            gauss_vec(vx) = water_vec{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = gauss_vec(vx);
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ex = nanmean(gauss_vec);
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.sd = nanstd(gauss_vec);
        end
    end


    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    ModelProcedureCell(mask_vec==0)=[];
    mask_vec(mask_vec==0)=[];
    
    
    final_water = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_water_sdex = toc(tstart);
    final_water_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_water{vx};
    end
    final_water_sdex = final_water_temp;
    save(fullfile([folder], 'final_water_sdex.mat'),'final_water_sdex','modeltime_final_water_sdex','-v7.3');
    delete(p);
end
%% Model metabolite data with metab cc prior knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)

for jobs = 1: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    
    mask_vec(mask_vec==0)=[];
    
    
    metab_cc = zeros(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    V = parallel.pool.Constant(mask_vec);
    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
        if V.Value(vx) >= 1
            [~, metab_cc(vx)] = fit_OspreyReferencing(D.Value{vx});
        end
    end
    metab_cc = metab_cc * data_vec{34}.txfrq*1e-6;
    modeltime_metab_cc = toc(tstart);
    save(fullfile([folder], 'metab_cc.mat'),'metab_cc','modeltime_metab_cc','-v7.3');
    delete(p);
    
    
    %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS.centerFreq = 4.68;
    BASIS = rmfield(BASIS,'specs');
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
     
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    ModelProcedureCell(mask_vec==0)=[];
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_cc(vx);
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = metab_cc(vx);
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = metab_cc(vx);
        end
    end
    
    
    
    
    final_metab_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_metab_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_metab_cc = toc(tstart);
    final_metab_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_metab_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_metab_cc{vx};
    end
    final_metab_cc = final_metab_cc_temp;
    save(fullfile([folder], 'final_metab_cc.mat'),'final_metab_cc','modeltime_final_metab_cc','-v7.3');
    delete(p);
end
%% Model metabolite data with metab no knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS);
    BASIS.centerFreq = 4.68;
    BASIS = rmfield(BASIS,'specs');
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
       
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    ModelProcedureCell(mask_vec==0)=[];
    mask_vec(mask_vec==0)=[];
    
    
    final_none = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_none(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_none = toc(tstart);
    final_none_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_none_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_none{vx};
    end
    final_none = final_none_temp;
    save(fullfile([folder], 'final_none.mat'),'final_none','modeltime_final_none','-v7.3');
    delete(p);
end

%% Model metabolite data with metab 5 cc prior knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                        'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                        'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    
    mask_vec(mask_vec==0)=[];
    
    load(fullfile([folder], 'metab_cc.mat'));
    
    
    
    %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS.centerFreq = 4.68;
    BASIS = rmfield(BASIS,'specs');
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_reduced_5_basis.json';
    
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    
    
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:5) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:5) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:5) = 1;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:5) = 0.5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:5) = 5;
    
    
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_cc(vx);
        end
    end
    
    ModelProcedureCell(mask_vec==0)=[];
    
    
    metab_5_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            metab_5_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_metab_5_cc = toc(tstart);
    save(fullfile([folder], 'metab_5_cc.mat'),'metab_5_cc','modeltime_metab_5_cc','-v7.3');
    delete(p);
    
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
    
    
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_5_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = metab_5_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = metab_5_cc{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end
    
    
    ModelProcedureCell(mask_vec==0)=[];
    
    final_metab_5_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_metab_5_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_metab_5_cc = toc(tstart);
    final_metab_5_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_metab_5_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_metab_5_cc{vx};
    end
    final_metab_5_cc = final_metab_5_cc_temp;
    save(fullfile([folder], 'final_metab_5_cc.mat'),'final_metab_5_cc','modeltime_final_metab_5_cc','-v7.3');
    delete(p);

    metab_5_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
    metab_5_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = metab_5_cc{vx};
    end
    metab_5_cc = metab_5_cc_temp;
    save(fullfile([folder], 'metab_5_cc.mat'),'metab_5_cc','modeltime_metab_5_cc','-v7.3');

end

%% Model metabolite data with metab mets cc prior knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1: length(matFiles)-4
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                        'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                        'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];
    
    mask_vec(mask_vec==0)=[];
    
    load(fullfile([folder], 'metab_cc.mat'));
    
    
    
    %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS = rmfield(BASIS,'specs');
    BASIS.centerFreq = 4.68;
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_reduced_mets_basis.json';
    
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    
    
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:18) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:18) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:18) = 1;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:18) = 0.5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:18) = 5;

    
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_cc(vx);
        end
    end
    
    ModelProcedureCell(mask_vec==0)=[];
    
    
    metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_metab_mets_cc = toc(tstart);
    save(fullfile([folder], 'metab_mets_cc.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    delete(p);
    
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
    
    
    ModelProcedureCell = cell(1,length(data_vec));
    
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end
    
    
    ModelProcedureCell(mask_vec==0)=[];
    
    final_metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_final_metab_mets_cc = toc(tstart);
    final_metab_mets_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_metab_mets_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_metab_mets_cc{vx};
    end
    final_metab_mets_cc = final_metab_mets_cc_temp;
    save(fullfile([folder], 'final_metab_mets_cc.mat'),'final_metab_mets_cc','modeltime_final_metab_mets_cc','-v7.3');
    delete(p);

     metab_mets_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
    metab_mets_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = metab_mets_cc{vx};
    end
    metab_mets_cc = metab_mets_cc_temp;
    save(fullfile([folder], 'metab_mets_cc.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
end

%% Model metabolite data squentially
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1 : length(matFiles)
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    % load(fullfile([folder], 'scaleData.mat'));
    % mask = squeeze(mask(:,:,2));
    
    conv = cell(size(mask));
    
    [cx,cy] = find(mask(:,:,2)==2);
    
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(MRSCont.processed.A{1}, BASIS); 
   
    BASIS.centerFreq = 4.68;

        % Get data scale
    scaleData = max(real(MRSCont.processed.A{1}.specs(MRSCont.processed.A{1}.ppm > -2 & MRSCont.processed.A{1}.ppm < 10 ,:,:,2)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    save(fullfile([folder], 'scaleData.mat'),'scaleData','-v7.3');

     BASIS = rmfield(BASIS,'specs');

    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);

 
    
    
    tstart = tic;
    data =op_takeVoxel(MRSCont.processed.A{1},[cx cy,2]);
    [~, fwhm] = fit_OspreyReferencing(data);
    fwhm = fwhm * MRSCont.processed.A{1}.txfrq*1e-6;
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_reduced_5_basis.json';
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.gaussLB.init = fwhm;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:5) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:5) = 2.42;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:5) = 1;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:5) = 0.5;
    % ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:5) = 5;

    
    
    conv(cx,cy,2) = Osprey_gLCM(data,ModelProcedure,0,1,scaleData,0,0,BASIS);
    
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    
    
    ModelProcedure.Steps{1}.parametrizations.ph0.init = conv{cx,cy,2}.Model{1}.parsOut.ph0;
    ModelProcedure.Steps{1}.parametrizations.ph1.init = conv{cx,cy,2}.Model{1}.parsOut.ph1;
    ModelProcedure.Steps{1}.parametrizations.gaussLB.init = conv{cx,cy,2}.Model{1}.parsOut.gaussLB;
    
    conv(cx,cy,2) = Osprey_gLCM(data,ModelProcedure,0,1,scaleData,0,0,BASIS);
    
    
    gaussLBvec = conv{cx,cy,2}.Model{1}.parsOut.gaussLB;
    ph0vec = conv{cx,cy,2}.Model{1}.parsOut.ph0;
    ph1vec = conv{cx,cy,2}.Model{1}.parsOut.ph1;

   order = [2 3 1];
   for z = 1 : 3   
        for nVox = 1 : 252
            [x,y] = osp_spiral(nVox);
            x = x+cx
            y = y+cy
    
            if  x > 0 && y > 0 && x < size(mask,1) && y < size(mask,2) && mask(x,y,order(z))> 0 
                if ~((x == cx) && (y == cy) && (order(z)==2))
                    data =op_takeVoxel(MRSCont.processed.A{1},[x y,order(z)]);
                    
                    ModelProcedure.Steps{1}.parametrizations.ph0.init = median(ph0vec);
                    ModelProcedure.Steps{1}.parametrizations.ph1.init = median(ph1vec);
                    ModelProcedure.Steps{1}.parametrizations.gaussLB.init = median(gaussLBvec);
        
                    conv(x,y,order(z)) = Osprey_gLCM(data,ModelProcedure,0,1,scaleData,0,0,BASIS);
        
                    gaussLBvec = [gaussLBvec conv{x,y,order(z)}.Model{1, 1}.parsOut.gaussLB];
                    ph0vec = [ph0vec conv{x,y,order(z)}.Model{1, 1}.parsOut.ph0];
                    ph1vec = [ph1vec conv{x,y,order(z)}.Model{1, 1}.parsOut.ph1];
                end
            end
        end
   end
    modeltime_conventional= toc(tstart);
    save(fullfile([folder], 'conv.mat'),'conv','modeltime_conventional','-v7.3');
end

%% Model water data squentially
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1 : length(matFiles)
    load(matFiles{jobs}); 
    load(maskFiles{jobs});

    MRSCont.processed.w{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    % mask = squeeze(mask(:,:,2));
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

    p = parpool('Dell');
    parfor vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            [~, refFWHM] = osp_XReferencing(data_vec{vx},4.68,1,[0 9.36],0);
            data_vec{vx}.FWHM = refFWHM * data_vec{vx}.txfrq*1e-6;
            fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs'};
                data_vec{vx} = rmfield(data_vec{vx},fields);
        end
    end
    delete(p);
    
    
    % Create basis with matching dwelltime 
    load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    DataToModel = op_zeropad(data_vec{1},2);   
    BASIS = fit_resampleBasis(DataToModel, BASIS); 
    BASIS.centerFreq = 4.68;
    
    
    % Get data scale
    % scaleData = max(real(MRSCont.processed.A{1}.specs(MRSCont.processed.A{1}.ppm > -2 & MRSCont.processed.A{1}.ppm < 10 ,:,:,2)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    % save(fullfile([folder], 'scaleData.mat'),'scaleData','-v7.3');
load(fullfile([folder], 'scaleData.mat'));

    BASIS = rmfield(BASIS,'specs');
    ToKeep = [16,29,30,31];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 4;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);

    clear MRSCont
    
    % Model water data 
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_water_in_vivo.json';
    
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
    

    ModelProcedureCell(mask_vec==0)=[];

    water = cell(1,length(data_vec));
    
    final_water = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_water = toc(tstart);
    delete(p);
    water_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = water{vx};
    end
    water = water_temp;
    save(fullfile([folder], 'water.mat'),'water','modeltime_water','-v7.3');
end


%% Model metabolite data with metab mets cc prior knowledge
clear all
[matFiles,maskFiles] = pathDef(2);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)

for jobs = 1: length(matFiles)
    jobs
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    load(fullfile([folder], 'scaleData.mat'));
    % load(fullfile([folder], 'metab_mets_cc.mat'));
    
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                % temp_metab_mets_cc{ind} = metab_mets_cc{x,y,z};
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];

    % temp_metab_mets_cc(mask_vec==0)=[];
    % metab_mets_cc = temp_metab_mets_cc;
    mask_vec(mask_vec==0)=[];
    
    
    % metab_cc = zeros(1,length(data_vec));
    % p = parpool('Dell');
    % D = parallel.pool.Constant(data_vec);
    % V = parallel.pool.Constant(mask_vec);
    % 
    % tstart = tic;
    % parfor vx = 1:  length(data_vec)
    %     if V.Value(vx) >= 1
    %         [~, metab_cc(vx)] = fit_OspreyReferencing(D.Value{vx});
    %     end
    % end
    % metab_cc = metab_cc * data_vec{34}.txfrq*1e-6;
    % modeltime_metab_cc = toc(tstart);
    % save(fullfile([folder], 'metab_cc.mat'),'metab_cc','modeltime_metab_cc','-v7.3');
    % delete(p);

    load(fullfile([folder], 'metab_cc.mat'));
     %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS = rmfield(BASIS,'specs');
    BASIS.centerFreq = 4.68;
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    
    
    % ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_reduced_mets_basis.json';
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_Spline_invivo_Reg_Optim_MRSI.json';
    % 
    ModelProcedure = jsonToStruct(ModelProced);
    ModelProcedure.Steps = ModelProcedure.Steps(1:2);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end



    ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:18) = 2.42;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:18) = 2.42;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:18) = 1;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:18) = 0.5;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:18) = 5;



    ModelProcedureCell = cell(1,length(data_vec));

    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_cc(vx);
        end
    end

    ModelProcedureCell(mask_vec==0)=[];


    metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_metab_mets_cc = toc(tstart);
    save(fullfile([folder], 'metab_mets_cc_optim_lessMM.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    delete(p);
    
    %  load(fullfile([folder], 'metab_mets_cc_optim.mat'));
    % ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';
    % 
    % ModelProcedure = jsonToStruct(ModelProced);
    % if isstruct(ModelProcedure.Steps)
    %     ModelProcedureCell = cell(size(ModelProcedure.Steps));
    %     for ss = 1 : size(ModelProcedure.Steps,1)
    %         ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    %     end
    %     ModelProcedure.Steps = ModelProcedureCell;
    % end
    % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 100;
    % 
    % 
    % ModelProcedureCell = cell(1,length(data_vec));
    % 
    % m = 14;
    % 
    % for vx = 1:  length(data_vec)
    %     if mask_vec(vx) >= 1
    %         ModelProcedureCell{vx} = ModelProcedure;
    %         ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
    %         for oo = 1 : 20
    %             AIC(oo) = log(sum((real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo)).^2))) + 2 * m * real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.ed(oo))/length(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo));
    %         end
    %        [~,minAICidx] = min(AIC);
    %         OptimmalRegPar = metab_mets_cc{1, vx}.Model{1, 2}.Regularization.Lambda(minAICidx);
    % 
    % 
    %         % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=metab_mets_cc{vx}.Model{1,2}.Regularization.OptimalRegPar;  
    %         ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=OptimmalRegPar;  
    %         % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
    %         % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
    %     end
    % end
    % 
    % 
    % ModelProcedureCell(mask_vec==0)=[];
    % 
    % final_metab_mets_cc = cell(1,length(data_vec));
    % p = parpool('Dell');
    % D = parallel.pool.Constant(data_vec);
    % M = parallel.pool.Constant(ModelProcedureCell);
    % S = parallel.pool.Constant(scaleData);
    % 
    % tstart = tic;
    % parfor vx = 1:  length(data_vec)
    %         final_metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    % end
    % modeltime_final_metab_mets_cc = toc(tstart);
    % final_metab_mets_cc_temp = cell(size(mask));
    % for vx = 1:  length(data_vec)
    %    final_metab_mets_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_metab_mets_cc{vx};
    % end
    % final_metab_mets_cc = final_metab_mets_cc_temp;
    % save(fullfile([folder], 'final_metab_mets_cc_spline_optim_m14_ind_lor.mat'),'final_metab_mets_cc','modeltime_final_metab_mets_cc','-v7.3');
    % delete(p);
    % % 
    % %  metab_mets_cc_temp = cell(size(mask));
    % % for vx = 1:  length(data_vec)
    % % metab_mets_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = metab_mets_cc{vx};
    % % end
    % % metab_mets_cc = metab_mets_cc_temp;
    % % save(fullfile([folder], 'metab_mets_cc.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    
end


%% Model GABA metabolite using short TE Gaussian LW
clear all
[matFiles,maskFiles] = pathDef(3);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1:length(matFiles)
% for jobs = 6:6
    load(matFiles{jobs}); 
    
    load(maskFiles{jobs});
    string = matFiles{jobs};
    MRSCont.processed.diff1{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    convfolder = folder;
    convfolder = strrep(convfolder,'edited','conv');
    load(fullfile([convfolder], 'scaleData.mat'));
    load(fullfile([convfolder], 'final_metab_mets_cc_spline_optim_m14_ind_lor.mat'));
    shortTE = final_metab_mets_cc;

    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.diff1{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                tempshortTE{ind} = shortTE{x,y,z};
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end

    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];

    tempshortTE(mask_vec==0)=[];
    shortTE = tempshortTE;
    mask_vec(mask_vec==0)=[];

    load(which('BASIS_Philips_Edited_se_MRSI_PRESS_GABA68_wMM.mat'));
    % BASIS.centerFreq = 4.68;
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(MRSCont.processed.diff1{1}, BASIS);
    basisSetOff = BASIS;
    basisSetOff.fids = basisSetOff.fids(:,:,1);
    basisSetOff.specs = basisSetOff.specs(:,:,1);
    basisSetOff.sz = size(basisSetOff.specs);
    basisSetOff.ppm = basisSetOff.ppm';
    BASIS.fids = BASIS.fids(:,:,3);
    BASIS.specs = BASIS.specs(:,:,3);
    BASIS.sz = size(BASIS.specs);
    fitOpts.coMM3 = '1to1GABA';
    fitOpts.FWHMcoMM3 = 14;
    fitOpts.CrFactor = 1;
    [BASIS] = osp_addDiffMMPeaks(BASIS,basisSetOff,fitOpts,4.68);
       
    BASIS.centerFreq = 4.68;

    % ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_Spline_invivo_Reg_Optim_MRSI_GABA.json';
    % % 
    % ModelProcedure = jsonToStruct(ModelProced);
    % ModelProcedure.Steps = ModelProcedure.Steps(1:2);
    % if isstruct(ModelProcedure.Steps)
    %     ModelProcedureCell = cell(size(ModelProcedure.Steps));
    %     for ss = 1 : size(ModelProcedure.Steps,1)
    %         ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    %     end
    %     ModelProcedure.Steps = ModelProcedureCell;
    % end
    % 
    % ModelProcedureCell = cell(1,length(data_vec));
    % 
    % for vx = 1:  length(data_vec)
    %     if mask_vec(vx) >= 1
    %         ModelProcedureCell{vx} = ModelProcedure;
    %         ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = shortTE{vx}.Model{1, 1}.parsOut.gaussLB;
    %     end
    % end
    % 
    % ModelProcedureCell(mask_vec==0)=[];
    % 
    % 
    % metab_mets_cc = cell(1,length(data_vec));
    % p = parpool('Dell');
    % D = parallel.pool.Constant(data_vec);
    % M = parallel.pool.Constant(ModelProcedureCell);
    % S = parallel.pool.Constant(scaleData);
    % tstart = tic;
    % parfor vx = 1:  length(data_vec)
    %         metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    % end
    % modeltime_metab_mets_cc = toc(tstart);
    % save(fullfile([folder], 'metab_mets_cc_optim.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    % delete(p);


    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut_GABA.json';
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
    end
    % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;


    m = 12;
    load(fullfile([folder], 'metab_mets_cc_optim.mat'))
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = shortTE{vx}.Model{1, 1}.parsOut.gaussLB;
            for oo = 1 : 20
                AIC(oo) = log(sum((real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo)).^2))) + 2 * m * real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.ed(oo))/length(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo));
            end
           [~,minAICidx] = min(AIC);
            OptimmalRegPar = metab_mets_cc{1, vx}.Model{1, 2}.Regularization.Lambda(minAICidx);
            ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=OptimmalRegPar;  
        end
    end


    ModelProcedureCell(mask_vec==0)=[];

    conv = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    tstart = tic;
    parfor vx = 1:  length(data_vec)
            conv(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_conventional = toc(tstart);
    conv_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       conv_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = conv{vx};
    end
    conv = conv_temp;
    save(fullfile([folder], 'conv_shortTE_LW_spline_optim_m12.mat'),'conv','modeltime_conventional','-v7.3');
    delete(p);   
end
%%
% Model GABA metabolite using short TE Gaussian LW
clear all
[matFiles,maskFiles] = pathDef(3);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1:length(matFiles)
% for jobs = 6:6
    load(matFiles{jobs}); 
    
    load(maskFiles{jobs});
    string = matFiles{jobs};
    MRSCont.processed.diff1{1}.nucleus = {'1H'};
    folder = MRSCont.outputFolder;
    folder = strrep(folder,'/Volumes/T7Shield/working','D:\working');
    folder = strrep(folder,'/','\');
    convfolder = folder;
    convfolder = strrep(convfolder,'edited','conv');
    load(fullfile([convfolder], 'scaleData.mat'));
    load(fullfile([convfolder], 'final_metab_mets_cc_spline.mat'));
    shortTE = final_metab_mets_cc;

    % mask = squeeze(mask(:,:,2));
    
    conv = cell(size(mask));
    
    [cx,cy] = find(mask(:,:,2)==2);

    load(which('BASIS_Philips_Edited_se_MRSI_PRESS_GABA68_wMM.mat'));
    % BASIS.centerFreq = 4.68;
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(MRSCont.processed.diff1{1}, BASIS);
    basisSetOff = BASIS;
    basisSetOff.fids = basisSetOff.fids(:,:,1);
    basisSetOff.specs = basisSetOff.specs(:,:,1);
    basisSetOff.sz = size(basisSetOff.specs);
    basisSetOff.ppm = basisSetOff.ppm';
    BASIS.fids = BASIS.fids(:,:,3);
    BASIS.specs = BASIS.specs(:,:,3);
    BASIS.sz = size(BASIS.specs);
    fitOpts.coMM3 = '1to1GABA';
    fitOpts.FWHMcoMM3 = 14;
    fitOpts.CrFactor = 1;
    [BASIS] = osp_addDiffMMPeaks(BASIS,basisSetOff,fitOpts,4.68);
       
    BASIS.centerFreq = 4.68;

    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut_GABA.json';
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
    end
    % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;
 
    

    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.diff1{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                tempshortTE{ind} = shortTE{x,y,z};
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end

    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];

    tempshortTE(mask_vec==0)=[];
    shortTE = tempshortTE;
    mask_vec(mask_vec==0)=[];

    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
        end
    end
    

    ModelProcedureCell(mask_vec==0)=[];
    
    conv = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            conv(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_conventional = toc(tstart);
    conv_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       conv_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = conv{vx};
    end
    conv = conv_temp;
    save(fullfile([folder], 'conv_no_LW_spline.mat'),'conv','modeltime_conventional','-v7.3');
    delete(p);   
end

%% TSC study 
% Model metabolite data with metab mets cc prior knowledge
clear all
[matFiles,maskFiles] = pathDef(4);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)


for jobs = 1: length(matFiles)
% for jobs = 11:11
    load(matFiles{jobs}); 
    load(maskFiles{jobs});

    MRSCont.processed.w{1}.nucleus = {'1H'};
    [~,container_name,~]=fileparts(matFiles{jobs});
    folder = fullfile('D:\working\MRSI\derivatives-TSC\conv',container_name(1:4))
    mkdir(folder);
    % mask = squeeze(mask(:,:,2));
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

    p = parpool('Dell');
    parfor vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            [~, refFWHM] = osp_XReferencing(data_vec{vx},4.68,1,[0 9.36],0);
            data_vec{vx}.FWHM = refFWHM * data_vec{vx}.txfrq*1e-6;
            fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs'};
                data_vec{vx} = rmfield(data_vec{vx},fields);
        end
    end
    delete(p);
    
    
    % Create basis with matching dwelltime 
    load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    DataToModel = op_zeropad(data_vec{1},2);   
    BASIS = fit_resampleBasis(DataToModel, BASIS); 
    BASIS.centerFreq = 4.68;
    
    
    % Get data scale
    scaleData = max(real(MRSCont.processed.A{1}.specs(MRSCont.processed.A{1}.ppm > -2 & MRSCont.processed.A{1}.ppm < 10 ,:,:,2)),[],'all') / max(max(max(real(BASIS.specs(BASIS.ppm > -2 & BASIS.ppm < 10 ,:)))));
    save(fullfile([folder], 'scaleData.mat'),'scaleData','-v7.3');
    % load(fullfile([folder], 'scaleData.mat'));

    BASIS = rmfield(BASIS,'specs');
    ToKeep = [16,29,30,31];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 4;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);

    clear MRSCont
    
    % Model water data 
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_water_in_vivo.json';
    
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
    

    ModelProcedureCell(mask_vec==0)=[];

    water = cell(1,length(data_vec));
    
    final_water = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            water(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_water = toc(tstart);
    delete(p);
    water_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       water_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = water{vx};
    end
    water = water_temp;
    save(fullfile([folder], 'water.mat'),'water','modeltime_water','-v7.3');
end
% Done with water model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for jobs = 1:length(matFiles)
% for jobs = 11:11
    load(matFiles{jobs}); 
    load(maskFiles{jobs});
    MRSCont.processed.A{1}.nucleus = {'1H'};
    [~,container_name,~]=fileparts(matFiles{jobs});
    folder = fullfile('D:\working\MRSI\derivatives-TSC\conv',container_name(1:4))
    mkdir(folder);
    load(fullfile([folder], 'scaleData.mat'));
    % load(fullfile([folder], 'metab_mets_cc.mat'));
    
    % mask = squeeze(mask(:,:,2));
    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.A{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                % temp_metab_mets_cc{ind} = metab_mets_cc{x,y,z};
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end
    
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];

    % temp_metab_mets_cc(mask_vec==0)=[];
    % metab_mets_cc = temp_metab_mets_cc;
    mask_vec(mask_vec==0)=[];
    
    
    metab_cc = zeros(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    V = parallel.pool.Constant(mask_vec);

    tstart = tic;
    parfor vx = 1:  length(data_vec)
        if V.Value(vx) >= 1
            [~, metab_cc(vx)] = fit_OspreyReferencing(D.Value{vx});
        end
    end
    metab_cc = metab_cc * data_vec{34}.txfrq*1e-6;
    modeltime_metab_cc = toc(tstart);
    save(fullfile([folder], 'metab_cc.mat'),'metab_cc','modeltime_metab_cc','-v7.3');
    delete(p);

    load(fullfile([folder], 'metab_cc.mat'));
     %Create basis with matching dwelltime 
     load(which('BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA20_wMM_expMM.mat'));
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard 
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    BASIS = rmfield(BASIS,'specs');
    BASIS.centerFreq = 4.68;
    ToKeep = [2,3,6,7,10,11,12,13,14,16,17,18,19,20,21,22,23,25,27];
    BASIS.fids = BASIS.fids(:,ToKeep);
    BASIS.name = BASIS.name(ToKeep);
    BASIS.nMets = 19;
    BASIS.nMM = 0;
    BASIS.sz = size(BASIS.fids);
    clear MRSCont
    
    
    % ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_reduced_mets_basis.json';
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_Spline_invivo_Reg_Optim_MRSI.json';
    % 
    ModelProcedure = jsonToStruct(ModelProced);
    ModelProcedure.Steps = ModelProcedure.Steps(1:2);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end



    ModelProcedure.Steps{1}.parametrizations.lorentzLB.init(1,1:18) = 2.42;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.ex(1,1:18) = 2.42;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.sd(1,1:18) = 1;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.lb(1,1:18) = 0.5;
    ModelProcedure.Steps{1}.parametrizations.lorentzLB.ub(1,1:18) = 5;



    ModelProcedureCell = cell(1,length(data_vec));

    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_cc(vx);
        end
    end

    ModelProcedureCell(mask_vec==0)=[];


    metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
            metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_metab_mets_cc = toc(tstart);
    save(fullfile([folder], 'metab_mets_cc_optim.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    delete(p);
    
     load(fullfile([folder], 'metab_mets_cc_optim.mat'));
    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut.json';
    % ModelProced = '/Volumes/Samsung/working/MRSI-model/model-procedures-expMM/1Step_Spline_invivo_Reg_Optim_MRSI.json';

    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 100;


    ModelProcedureCell = cell(1,length(data_vec));

    m = 14;

    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            for oo = 1 : 20
                AIC(oo) = log(sum((real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo)).^2))) + 2 * m * real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.ed(oo))/length(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo));
            end
           [~,minAICidx] = min(AIC);
            OptimmalRegPar = metab_mets_cc{1, vx}.Model{1, 2}.Regularization.Lambda(minAICidx);


            % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=metab_mets_cc{vx}.Model{1,2}.Regularization.OptimalRegPar;  
            ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=OptimmalRegPar;  
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.lb = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
            % ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.ub = metab_mets_cc{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end


    ModelProcedureCell(mask_vec==0)=[];

    final_metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    tstart = tic;
    parfor vx = 1:  length(data_vec)
            final_metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_final_metab_mets_cc = toc(tstart);
    final_metab_mets_cc_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       final_metab_mets_cc_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = final_metab_mets_cc{vx};
    end
    final_metab_mets_cc = final_metab_mets_cc_temp;
    save(fullfile([folder], 'final_metab_mets_cc_spline_optim_m14_ind_lor.mat'),'final_metab_mets_cc','modeltime_final_metab_mets_cc','-v7.3');
    delete(p);
    
end


%% Model GABA metabolite using short TE Gaussian LW
clear all
[matFiles,maskFiles] = pathDef(5);
%maskFiles=maskFiles(7)
%matFiles=matFiles(7)
% for jobs = 1:length(matFiles)
for jobs = 18:18
    load(matFiles{jobs}); 
    
    load(maskFiles{jobs});
    string = matFiles{jobs};
    MRSCont.processed.diff1{1}.nucleus = {'1H'};
    [~,container_name,~]=fileparts(matFiles{jobs});
    folder = fullfile('D:\working\MRSI\derivatives-TSC\edited',container_name(1:4))
    mkdir(folder);
    convfolder = folder;
    convfolder = strrep(convfolder,'edited','conv');
    load(fullfile([convfolder], 'scaleData.mat'));
    load(fullfile([convfolder], 'final_metab_mets_cc_spline_optim_m14_ind_lor.mat'));
    shortTE = final_metab_mets_cc;

    x_vec=[];
    y_vec=[];
    z_vec=[];
    % Unpack data into list
    ind = 1;
    for z = 1 : MRSCont.processed.A{1}.nZvoxels
        for x = 1 : MRSCont.processed.A{1}.nXvoxels
            for y = 1 : MRSCont.processed.A{1}.nYvoxels
                data_vec{ind}=op_takeVoxel(MRSCont.processed.diff1{1},[x y z]);
                fields = {'nXvoxels','nYvoxels','nZvoxels','seq','PatientPosition',...
                            'Manufacturer','OriginalFile','software','nii_mrs','specReg','watersupp',...
                            'refFWHM','refShift'};
                data_vec{ind} = rmfield(data_vec{ind},fields);
                mask_vec(ind) = mask(x,y,z);
                tempshortTE{ind} = shortTE{x,y,z};
                ind = ind + 1;
                x_vec=[ x_vec x];
                y_vec=[ y_vec y];
                z_vec=[ z_vec z];
            end
        end
    end

    if jobs == 18
        mask_vec(268:end) = 0;
    end
    data_vec(mask_vec==0)=[];
    x_vec(mask_vec==0)=[];
    y_vec(mask_vec==0)=[];
    z_vec(mask_vec==0)=[];

    tempshortTE(mask_vec==0)=[];
    shortTE = tempshortTE;
    mask_vec(mask_vec==0)=[];

    load(which('BASIS_Philips_Edited_se_MRSI_PRESS_GABA68_wMM.mat'));
    % BASIS.centerFreq = 4.68;
    BASIS = recalculateBasisSpecs(BASIS);                         % Add ppm axis and frequency domain data
    BASIS = fit_sortBasisSet(BASIS);                              % Sort according to Osprey standard
    BASIS = fit_resampleBasis(data_vec{1}, BASIS); 
    basisSetOff = BASIS;
    basisSetOff.fids = basisSetOff.fids(:,:,1);
    basisSetOff.specs = basisSetOff.specs(:,:,1);
    basisSetOff.sz = size(basisSetOff.specs);
    basisSetOff.ppm = basisSetOff.ppm';
    BASIS.fids = BASIS.fids(:,:,3);
    BASIS.specs = BASIS.specs(:,:,3);
    BASIS.sz = size(BASIS.specs);
    fitOpts.coMM3 = '1to1GABA';
    fitOpts.FWHMcoMM3 = 14;
    fitOpts.CrFactor = 1;
    [BASIS] = osp_addDiffMMPeaks(BASIS,basisSetOff,fitOpts,4.68);
       
    BASIS.centerFreq = 4.68;

    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_Spline_invivo_Reg_Optim_MRSI_GABA.json';

    ModelProcedure = jsonToStruct(ModelProced);
    ModelProcedure.Steps = ModelProcedure.Steps(1:2);
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end

    ModelProcedureCell = cell(1,length(data_vec));

    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = shortTE{vx}.Model{1, 1}.parsOut.gaussLB;
        end
    end

    ModelProcedureCell(mask_vec==0)=[];


    metab_mets_cc = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);
    tstart = tic;
    parfor vx = 1:  length(data_vec)
        vx
            metab_mets_cc(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,0,S.Value,0,0,BASIS);
    end
    modeltime_metab_mets_cc = toc(tstart);
    save(fullfile([folder], 'metab_mets_cc_optim.mat'),'metab_mets_cc','modeltime_metab_mets_cc','-v7.3');
    delete(p);


    ModelProced = 'C:\Users\Superuser\Documents\GitHub\MRSI-model\model-procedures\1Step_final_in_vivo_cut_GABA.json';
    ModelProcedure = jsonToStruct(ModelProced);
    if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
    end
    % ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 15;


    m = 12;
    load(fullfile([folder], 'metab_mets_cc_optim.mat'))
    for vx = 1:  length(data_vec)
        if mask_vec(vx) >= 1
            ModelProcedureCell{vx} = ModelProcedure;
            ModelProcedureCell{vx}.Steps{1}.parametrizations.gaussLB.init = shortTE{vx}.Model{1, 1}.parsOut.gaussLB;
            for oo = 1 : 20
                AIC(oo) = log(sum((real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo)).^2))) + 2 * m * real(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.ed(oo))/length(metab_mets_cc{1, vx}.Model{1, 2}.Regularization.residual(:,oo));
            end
           [~,minAICidx] = min(AIC);
            OptimmalRegPar = metab_mets_cc{1, vx}.Model{1, 2}.Regularization.Lambda(minAICidx);
            ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar=OptimmalRegPar;  
        end
    end


    ModelProcedureCell(mask_vec==0)=[];

    conv = cell(1,length(data_vec));
    p = parpool('Dell');
    D = parallel.pool.Constant(data_vec);
    M = parallel.pool.Constant(ModelProcedureCell);
    S = parallel.pool.Constant(scaleData);

    tstart = tic;
    parfor vx = 1:  length(data_vec)
        vx
            conv(vx) = Osprey_gLCM(D.Value(vx),M.Value{vx},0,1,S.Value,0,0,BASIS);
    end
    modeltime_conventional = toc(tstart);
    conv_temp = cell(size(mask));
    for vx = 1:  length(data_vec)
       conv_temp{x_vec(vx),y_vec(vx),z_vec(vx)} = conv{vx};
    end
    conv = conv_temp;
    save(fullfile([folder], 'conv_shortTE_LW_spline_optim_m12.mat'),'conv','modeltime_conventional','-v7.3');
    delete(p);   
end
