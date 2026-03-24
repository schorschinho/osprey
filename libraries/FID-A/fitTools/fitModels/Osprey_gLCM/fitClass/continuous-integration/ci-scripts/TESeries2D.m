classdef TESeries2D < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',0,'plot',0);
        baseline = {'none','spline'};
        optimSignalPart = {'R'};
    end
    
    methods (Test)

        function test_2D_TESeries_sI(testCase,baseline,optimSignalPart)            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/2D/TE-SeriesNoBsNoMMIns/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/TESeries2D/TE_Series_1Step_Spline_sI.json');
           
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
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.5;
                case 'poly'
                    ModelProcedure.Steps{1}.fit_opts.baseline.order = 4;
                case 'reg_spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.type = 'spline';
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.0667;
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegFun = 'SecondDerivative';
                    ModelProcedure.Steps{1}.parametrizations.baseAmpl.RegPar = 500;           
            end

            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
            
            % Verify amplitude estimates
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl(1);
            ampl_gt = 0.15;

            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_sI_ampl_reg_spline_' optimSignalPart '.mat']);
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