classdef BetweenStepParser < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',0,'plot',0);
    end
    
    methods (Test)
        
        function test_1D_3Steps(testCase)   
            
            % Fit with parameterized MMs
            % Find data and model procedure
            container = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/jobSDAT.mat');
            load(container);
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
            
            % Set inital Gauss LB value
            MRSCont.processed.metab{1}.FWHM = MRSCont.QM.FWHM.metab(1);

            % Run gLCM
            test = Osprey_gLCM(MRSCont.processed.metab{1},ModelProcedure,0,1);
            
            % Verify how data is parsed through the 3 steps
            verifyEqual(testCase, test{1}.step, 3, 'Did we run 3 steps?')   
            % Validate step 1
            verifyEqual(testCase, test{1}.Options{1}.optimFreqFitRange, [1.8; 4], 'Did we use the right fit range?') 
            verifyEqual(testCase,test{1}.Options{1}.parametrizations.gaussLB.init, MRSCont.processed.metab{1}.FWHM , 'Did we use inital Gauss?')
            verifyEqual(testCase,test{1}.Options{1}.FunctionTolerance,1e-3, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{1}.StepTolerance,1e-3, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{1}.OptimalityTolerance,1e-3, 'Did we correct tolerances?')
            % Validate step 2
            verifyEqual(testCase, test{1}.Options{2}.optimFreqFitRange, [0.2; 4], 'Did we use the right fit range?') 
            params = {'ph0','ph1','gaussLB'};
            for pp = 1 : length(params)
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).init, test{1}.Model{1}.parsOut.(params{pp}) , ['Did we use inital' params{pp} ' from step 1?'])
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).ub, test{1}.Model{1}.parsOut.(params{pp}) , ['Did we update ub' params{pp} ' from step 1?'])
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).lb, test{1}.Model{1}.parsOut.(params{pp}) , ['Did we update lb' params{pp} ' from step 1?'])
            end
            params = {'freqShift','lorentzLB'};
            for pp = 1 : length(params)
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).init(1:18), test{1}.Model{1}.parsOut.(params{pp}) , ['Did we use inital' params{pp} ' from step 1?'])
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).ub(1:18), test{1}.Model{1}.parsOut.(params{pp}) , ['Did we update ub' params{pp} ' from step 1?'])
                verifyEqual(testCase,test{1}.Options{2}.parametrizations.(params{pp}).lb(1:18), test{1}.Model{1}.parsOut.(params{pp}) , ['Did we update lb' params{pp} ' from step 1?'])
            end
            verifyEqual(testCase,test{1}.Options{2}.FunctionTolerance,1e-3, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{2}.StepTolerance,1e-3, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{2}.OptimalityTolerance,1e-3, 'Did we correct tolerances?')
            % Validate step 2
            verifyEqual(testCase,test{1}.Options{3}.parametrizations.baseAmpl.RegPar, test{1}.Model{2}.Regularization.OptimalRegPar , ['Did we use the correct regularization parameter?'])
            verifyEqual(testCase,test{1}.Options{3}.FunctionTolerance,1e-10, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{3}.StepTolerance,1e-10, 'Did we use correct tolerances?')
            verifyEqual(testCase,test{1}.Options{3}.OptimalityTolerance,1e-10, 'Did we correct tolerances?')
            
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