classdef BaselineModels1D < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',0,'plot',0);
        baseline = {'none','spline','poly','reg_spline'};
        optimSignalPart = {'R','RI'};
    end
    
    methods (Test)

        function test_1D_sI_singlestep_baseline_models(testCase,baseline,optimSignalPart)            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/sINoBS/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_sI.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline;
            switch baseline
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.15;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end

            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl;
            ampl_gt = 0.15;

            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_sI.json',['ci-prior-iteration/BaselineModels1D/1Step_sI_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 0.1)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.01)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end

        function test_1D_invivo_noBS_noMM_singlestep_baseline_models(testCase,baseline,optimSignalPart)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/NoBSNoMM/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_invivo.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline;
            switch baseline
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.0667;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end
            ModelProcedure.Steps{1}.basisset.include = ModelProcedure.Steps{1}.basisset.include(1:18);
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39];

            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_noMM_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 2.5)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.001)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end

        function test_1D_invivo_noBS_singlestep_baseline_models(testCase,baseline,optimSignalPart)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/NoBS/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_invivo.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline;
            switch baseline
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.0667;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;
            ModelProcedure.basisset.mmdef = [];

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39,...
                       1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];


            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_noBS_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 10)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.002)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end

        function test_1D_invivo_singlestep_baseline_models(testCase,baseline,optimSignalPart)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_invivo.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline;
            switch baseline
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.0667;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;
            ModelProcedure.basisset.mmdef = [];

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39,...
                       1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];

            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            if ~strcmp(baseline,'none') % Makes no sense to compare this
                verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 20)
            end
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.01)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end

        function test_1D_invivo_ExpMM_singlestep_baseline_models(testCase,baseline,optimSignalPart)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/ExpMM/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_invivo.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline;
            switch baseline
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.0667;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end
            ModelProcedure.Steps{1}.basisset.include = ModelProcedure.Steps{1}.basisset.include(1:19);
            ModelProcedure.Steps{1}.basisset.include{19} = 'MM_PRESS_PCC';
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;
            ModelProcedure.basisset.mmdef = [];
            
            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39,...
                       7.2];

            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_invivo.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 3.4)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.001)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end

        function test_1D_invivo_optimized_baseline_models(testCase,optimSignalPart)   
            
            % Fit with parameterized MMs
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/BaselineModels1D/1Step_Spline_invivo_Reg_Optim_Full.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end
            ModelProcedure.basisset.mmdef = [];

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{3}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39,...
                       1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];


            % Save outputs for comparison
            if testCase.debug.update_priors
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_Spline_invivo_Reg_Optim_Full.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_BSandMMs_optim_spline_ampl_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
            else
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_Spline_invivo_Reg_Optim_Full.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_BSandMMs_optim_spline_ampl_' optimSignalPart '.mat']);
                        load(PriorFile);
            end

            % Verify with gt values
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 10)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.002)

            if testCase.debug.plot
                test{1}.plotFit(3)
            end

            % Fit with exp MMs
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/ExpMM/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');

            % Parse parameters
            ModelProcedure.Steps{2}.basisset.include = ModelProcedure.Steps{2}.basisset.include(1:19);
            ModelProcedure.Steps{2}.basisset.include{19} = 'MM_PRESS_PCC';
            ModelProcedure.Steps{3}.basisset.include = ModelProcedure.Steps{3}.basisset.include(1:19);
            ModelProcedure.Steps{3}.basisset.include{19} = 'MM_PRESS_PCC';
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;
            ModelProcedure.basisset.mmdef = [];

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{3}.parsOut.metAmpl;
            ampl_gt = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39,...
                       7.2];

            % Save outputs for comparison
            if testCase.debug.update_priors
                       PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_Spline_invivo_Reg_Optim_Full.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_optim_spline_ampl_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
            else
                        PriorFile = strrep(model,'ci-model-procedures/BaselineModels1D/1Step_Spline_invivo_Reg_Optim_Full.json',['ci-prior-iteration/BaselineModels1D/1Step_invivo_ExpMM_optim_spline_ampl_' optimSignalPart '.mat']);
                        load(PriorFile);
            end

            % Verify with gt values
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 10)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.002)

            if testCase.debug.plot
                test{1}.plotFit(3)
            end
        end

    end
    
    methods (TestClassSetup)       
        function classSetup(testCase)
            % Set up shared state for all tests. 
            import matlab.unittest.constraints.IsEqualTo
            import matlab.automation.diagnostics.DisplayDiagnostic
        
            % Find folder with data 
            DataFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/'); 
            addpath(genpath(DataFolder));                                               % Add data folder to path
            % Find folder with prior results 
            PriorFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'ci-prior-iteration/'); 
            addpath(genpath(PriorFolder));                                               % Add prior results folder to path
        end     
    end

    methods (TestMethodTeardown)
        function deleteData(testCase)

        end
    end
    
    % methods (Test, ParameterCombination = 'sequential')
    %     function testNumel(testCase,level,side)
    %         import matlab.unittest.constraints.HasElementCount
    %         testCase.verifyThat(sierpinski(level), ...
    %             HasElementCount(side^2))
    %     end
    % end 
end