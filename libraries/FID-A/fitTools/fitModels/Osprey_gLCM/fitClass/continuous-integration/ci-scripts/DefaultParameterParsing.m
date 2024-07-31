classdef DefaultParameterParsing < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',0,'plot',0);
        parameter = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'};
    end
    
    methods (Test)

        

        function test_1D_default_parameter_parsing(testCase,parameter)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/NIfTIMRS/sim-001_ses-001_MRS_3T_30_TE_A.nii.gz');
            model = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/ParameterParsing/1Step_invivo_default.json');
           
            % Load model json
            ModelProcedure = jsonToStruct(model);
            if isstruct(ModelProcedure.Steps)
                ModelProcedureCell = cell(size(ModelProcedure.Steps));
                for ss = 1 : size(ModelProcedure.Steps,1)
                    ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
                end
                ModelProcedure.Steps = ModelProcedureCell;
            end
            
            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);

            % Set default initals
            % Initialize phi0 as constant with value 0
            defaults.ph0.fun    = 'free';
            defaults.ph0.gradfun = 'free';
            defaults.ph0.lb      = -pi/4;
            defaults.ph0.ub      = pi/4;
            defaults.ph0.init    = 0;
            defaults.ph0.ex      = 0;
            defaults.ph0.sd      = Inf;
            defaults.ph0.RegFun  = '';
            defaults.ph0.gr      = [];

            % Initialize phi1 as constant with value 0
            defaults.ph1.fun     = 'free';
            defaults.ph1.gradfun = 'free';
            defaults.ph1.lb      = -pi/270;
            defaults.ph1.ub      = pi/270;
            defaults.ph1.init    = 0;
            defaults.ph1.ex      = 0;
            defaults.ph1.sd      = Inf;
            defaults.ph1.RegFun  = '';
            defaults.ph1.gr      = [];
    
            % Initialize Gaussian LB as constant with value [0.04 *
            % hz/ppm]
            defaults.gaussLB.fun     = 'free';
            defaults.gaussLB.gradfun = 'free';
            defaults.gaussLB.lb      = 0;
            defaults.gaussLB.ub      = 70;
            defaults.gaussLB.init    = 0.04 * test{1}.Data.txfrq*1e-6;
            defaults.gaussLB.ex      = 0.04 * test{1}.Data.txfrq*1e-6;
            defaults.gaussLB.sd      = Inf;
            defaults.gaussLB.RegFun  = '';
            defaults.gaussLB.gr      = [];
    
            % Initialize Lorentzian LB as constant with value
            defaults.lorentzLB.fun     = 'free';
            defaults.lorentzLB.gradfun = 'free';
            defaults.lorentzLB.lb      = 0.5*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.lorentzLB.ub      = 7*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.lorentzLB.init    = 2.75*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.lorentzLB.ex      = 2.75*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.lorentzLB.sd      = 2.75*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.lorentzLB.RegFun  = '';
            defaults.lorentzLB.gr      = [];
    
            % Initialize frequency shifts as constant with value 0 Hz
            defaults.freqShift.fun     = 'free';
            defaults.freqShift.gradfun = 'free';
            defaults.freqShift.lb      = -10*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.freqShift.ub      = 10*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.freqShift.init    = zeros(1, sum(test{1}.BasisSets.includeInFit));
            defaults.freqShift.ex      = zeros(1, sum(test{1}.BasisSets.includeInFit));
            defaults.freqShift.sd      = 0.04 * test{1}.Data.txfrq*1e-6 *ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.freqShift.RegFun  = '';
            defaults.freqShift.gr      = [];
    
            % Initialize metabolite amplitudes as free with value 0
            defaults.metAmpl.fun     = 'free';
            defaults.metAmpl.gradfun = 'free';
            defaults.metAmpl.lb      = zeros(1, sum(test{1}.BasisSets.includeInFit));
            defaults.metAmpl.ub      = Inf*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.metAmpl.init    = zeros(1, sum(test{1}.BasisSets.includeInFit));
            defaults.metAmpl.ex      = zeros(1, sum(test{1}.BasisSets.includeInFit));
            defaults.metAmpl.sd      = Inf*ones(1, sum(test{1}.BasisSets.includeInFit));
            defaults.metAmpl.RegFun  = '';
            defaults.metAmpl.gr      = [];
    
            % Initialize baseline amplitudes as free with value 0
            defaults.baseAmpl.fun     = 'free';
            defaults.baseAmpl.gradfun = 'free';
            defaults.baseAmpl.lb      = -Inf*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.ub      = Inf*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.init    = zeros(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.ex      = zeros(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.sd      = Inf*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.RegFun  = '';
            defaults.baseAmpl.gr      = [];

            field_names = fieldnames(defaults.ph0);
            for ff = 1 : length(field_names)
                verifyEqual(testCase,test{1}.Options{1}.parametrizations.(parameter).(field_names{ff}), defaults.(parameter).(field_names{ff}) , ['Did set the default ' field_names{ff} ' from ' parameter ' correctly?'])
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