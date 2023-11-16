classdef MultiTransient2D < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',1,'plot',0);
        baseline = {'none','spline','poly','reg_spline'};
        baseline_red = {'spline'};
        optimSignalPart = {'R'};
        parameter = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'};
        type = {'fixed','free'};
    end
    
    methods (Test)

        function test_2D_sI_singlestep_baseline_models(testCase,baseline,optimSignalPart)            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/2D/sINoBSSeparateTransients/NIfTIMRS/sim-001_ses-001_MRS_3T_30_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_sI.json');
           
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

      
        function test_2D_NAAmICr_noBS_singlestep_baseline_models(testCase,baseline,optimSignalPart)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/2D/NoBSNoMMNAACrInsSeparateTransients/NIfTIMRS/sim-001_ses-001_MRS_3T_30_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json');
           
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
            ampl = test{1, 1}.Model{1, 1}.parsOut.metAmpl(1,:);
            ampl_gt = [3.13,7.06,8.70];


            % Save outputs for comparison
            if testCase.debug.update_priors
               switch baseline
                    case 'none'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_none_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'spline'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                    case 'poly'
                       PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_poly_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version'); 
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_reg_spline_' optimSignalPart '.mat']);
                       ampl_prior_version = round(ampl,4);
                       save(PriorFile,'ampl_prior_version');
                       baseline = 'spline';
                    end
            else
                switch baseline
                    case 'none'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_none_' optimSignalPart '.mat']);
                        load(PriorFile);
                    case 'spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_spline_' optimSignalPart '.mat']);
                        load(PriorFile); 
                    case 'poly'
                         PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_poly_' optimSignalPart '.mat']);
                         load(PriorFile);
                    case 'reg_spline'
                        PriorFile = strrep(model,'ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json',['ci-prior-iteration/MultiTransientModels2D/2D_1Step_invivo_noBS_ampl_reg_spline_' optimSignalPart '.mat']);
                        load(PriorFile);   
                        baseline = 'spline';
                end
            end

            % Verify with gt values
            verifyEqual(testCase, test{1}.Options{1}.baseline.type, baseline, 'Does the baseline model match?')
            verifyEqual(testCase, ampl, ampl_gt, 'Does the estimate match the groundtruth?', 'RelTol', 10)
            verifyEqual(testCase, ampl, ampl_prior_version, 'Does the estimate match the prior version of the gLCM?', 'RelTol', 0.02)

            if testCase.debug.plot
                test{1}.plotFit(1)
            end
        end  

        function test_2D_parametrization_multi_transient(testCase,baseline_red,optimSignalPart,parameter,type)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/2D/NoBSNoMMNAACrInsSeparateTransients/NIfTIMRS/sim-001_ses-001_MRS_3T_30_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/MultiTransientModels2D/2D_1Step_NoBS_NAAmICr.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end

            % Load dynamic model json
            DynamicModelJson = ConvertRelativePath(ModelProcedure.Steps{1}.extra.DynamicModelJson);
            DynamicModel = jsonToStruct(DynamicModelJson);
            ModelProcedure.Steps{1}.extra.DynamicModelJson = DynamicModel;

            % Parse parameters
            ModelProcedure.Steps{1}.fit_opts.baseline.type = baseline_red;
            switch baseline_red
                case 'spline'
                    ModelProcedure.Steps{1}.fit_opts.baseline.dkntmn = 0.4;    
            end
            ModelProcedure.Steps{1, 1}.fit_opts.optimSignalPart = optimSignalPart;
            ModelProcedure.basisset.mmdef = [];

            switch parameter
                case {'ph0', 'ph1', 'gaussLB'}
                    nPars = 1;
                case {'lorentzLB','freqShift','metAmpl'}
                    nPars = length(ModelProcedure.Steps{1}.basisset.include);
                case 'baseAmpl'
                    nPars = 11;
            end

                       
            % Set default initals and parametrizations
            % Initialize phi0 as constant with value 0
            defaults.ph0.fun    = 'free';
            defaults.ph0.gradfun = 'free';
            defaults.ph0.lb      = -pi/4;
            defaults.ph0.ub      = pi/4;
            defaults.ph0.init    = 0;
            defaults.ph0.ex      = 0;
            defaults.ph0.sd      = Inf;
            defaults.ph0.RegFun  = '';
            defaults.ph0.type    = 'fixed';

            % Initialize phi1 as constant with value 0
            defaults.ph1.fun     = 'free';
            defaults.ph1.gradfun = 'free';
            defaults.ph1.lb      = -pi/270;
            defaults.ph1.ub      = pi/270;
            defaults.ph1.init    = 0;
            defaults.ph1.ex      = 0;
            defaults.ph1.sd      = Inf;
            defaults.ph1.RegFun  = '';
            defaults.ph1.type    = 'fixed';
    
            % Initialize Gaussian LB as constant with value [0.04 *
            % hz/ppm]
            defaults.gaussLB.fun     = 'free';
            defaults.gaussLB.gradfun = 'free';
            defaults.gaussLB.lb      = 0;
            defaults.gaussLB.ub      = 70;
            defaults.gaussLB.init    = 9.70;
            defaults.gaussLB.ex      = 9.70;
            defaults.gaussLB.sd      = Inf;
            defaults.gaussLB.RegFun  = '';
            defaults.gaussLB.type    = 'fixed';
    
            % Initialize Lorentzian LB as constant with value
            defaults.lorentzLB.fun     = 'free';
            defaults.lorentzLB.gradfun = 'free';
            defaults.lorentzLB.lb      = 0.5*ones(1, nPars);
            defaults.lorentzLB.ub      = 7*ones(1, nPars);
            defaults.lorentzLB.init    = 3.40*ones(1, nPars);
            defaults.lorentzLB.ex      = 2.75*ones(1, nPars);
            defaults.lorentzLB.sd      = 2.75*ones(1, nPars);
            defaults.lorentzLB.RegFun  = '';
            defaults.lorentzLB.type    = 'fixed';
    
            % Initialize frequency shifts as constant with value 0 Hz
            defaults.freqShift.fun     = 'free';
            defaults.freqShift.gradfun = 'free';
            defaults.freqShift.lb      = -10*ones(1, nPars);
            defaults.freqShift.ub      = 10*ones(1, nPars);
            defaults.freqShift.init    = zeros(1, nPars);
            defaults.freqShift.ex      = zeros(1, nPars);
            defaults.freqShift.sd      = 3*ones(1, nPars);
            defaults.freqShift.RegFun  = '';
            defaults.freqShift.type    = 'fixed';
    
            % Initialize metabolite amplitudes as free with value 0
            defaults.metAmpl.fun     = 'free';
            defaults.metAmpl.gradfun = 'free';
            defaults.metAmpl.lb      = zeros(1, nPars);
            defaults.metAmpl.ub      = Inf*ones(1, nPars);
            defaults.metAmpl.init    = [3.13,7.06,8.70];
            defaults.metAmpl.ex      = zeros(1, nPars);
            defaults.metAmpl.sd      = Inf*ones(1, nPars);
            defaults.metAmpl.RegFun  = '';
            defaults.metAmpl.type    = 'fixed';
    
            % Initialize baseline amplitudes as free with value 0
            defaults.baseAmpl.fun     = 'free';
            defaults.baseAmpl.gradfun = 'free';
            defaults.baseAmpl.lb      = -Inf*ones(1, nPars);
            defaults.baseAmpl.ub      = Inf*ones(1, nPars);
            defaults.baseAmpl.init    = zeros(1, nPars);
            defaults.baseAmpl.ex      = zeros(1, nPars);
            defaults.baseAmpl.sd      = Inf*ones(1, nPars);
            defaults.baseAmpl.RegFun  = '';
            defaults.baseAmpl.type    = 'fixed';
        
            % Update other fields
            ModelProcedure.Steps{1}.extra.DynamicModelJson.parameters.(parameter).type = type;
            
            % Overwrite with expected values according to parametrization
             defaults.(parameter).type = type;
             fields = {'lb','ub','init','ex','sd'};
             for ff = 1 : length(fields)
                 if strcmp(type,'free')
                    defaults.(parameter).(fields{ff}) = repmat(defaults.(parameter).(fields{ff}), [16 1]);
                 end
             end

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);

            field_names = fieldnames(defaults.ph0);
            for ff = 1 : length(field_names)
                verifyEqual(testCase,test{1}.Options{1}.parametrizations.(parameter).(field_names{ff}), defaults.(parameter).(field_names{ff}) , ['Did the reparametrization of ' field_names{ff} ' from ' parameter ' work correctly?'])
            end
                      
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