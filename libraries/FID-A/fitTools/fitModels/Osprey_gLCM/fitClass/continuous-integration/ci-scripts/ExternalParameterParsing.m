classdef ExternalParameterParsing < matlab.unittest.TestCase   
    properties (TestParameter)
        debug = struct('update_priors',0,'plot',0);
        parameter = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'};
        parameter_per_basis = {'lorentzLB_indiv','freqShift_indiv','metAmpl_indiv','baseAmpl_indiv'};
        fields = {'lb','ub','init','ex','sd'}
    end
    
    methods (Test)

        

        function test_1D_external_parameter_parsing(testCase,parameter)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which(['libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/ParameterParsing/1Step_invivo_external_' parameter '.json']);
           
            
            
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
            
            nPars = length(ModelProcedure.Steps{1}.basisset.include);

            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);

            % Set default initals
            % Initialize phi0 as constant with value 0
            defaults.ph0.fun    = 'free';
            defaults.ph0.gradfun = 'free';
            defaults.ph0.lb      = -1;
            defaults.ph0.ub      = 1;
            defaults.ph0.init    = 0.1;
            defaults.ph0.ex      = 0.01;
            defaults.ph0.sd      = 100;
            defaults.ph0.RegFun  = '';

            % Initialize phi1 as constant with value 0
            defaults.ph1.fun     = 'free';
            defaults.ph1.gradfun = 'free';
            defaults.ph1.lb      = -1;
            defaults.ph1.ub      = 1;
            defaults.ph1.init    = 0.1;
            defaults.ph1.ex      = 0.01;
            defaults.ph1.sd      = 100;
            defaults.ph1.RegFun  = '';
    
            % Initialize Gaussian LB as constant with value [0.04 *
            % hz/ppm]
            defaults.gaussLB.fun     = 'free';
            defaults.gaussLB.gradfun = 'free';
            defaults.gaussLB.lb      = 1;
            defaults.gaussLB.ub      = 35;
            defaults.gaussLB.init    = 3;
            defaults.gaussLB.ex      = 3;
            defaults.gaussLB.sd      = 100;
            defaults.gaussLB.RegFun  = '';

            % Initialize Lorentzian LB as constant with value
            defaults.lorentzLB.fun     = 'free';
            defaults.lorentzLB.gradfun = 'free';
            defaults.lorentzLB.lb      = ones(1, nPars);
            defaults.lorentzLB.ub      = 8*ones(1, nPars);
            defaults.lorentzLB.init    = 3.4*ones(1, nPars);
            defaults.lorentzLB.ex      = 3.4*ones(1, nPars);
            defaults.lorentzLB.sd      = 1000*ones(1, nPars);
            defaults.lorentzLB.RegFun  = '';
    
            % Initialize frequency shifts as constant with value 0 Hz
            defaults.freqShift.fun     = 'free';
            defaults.freqShift.gradfun = 'free';
            defaults.freqShift.lb      = -5*ones(1, nPars);
            defaults.freqShift.ub      = 5*ones(1, nPars);
            defaults.freqShift.init    = ones(1, nPars);
            defaults.freqShift.ex      = 2*ones(1, nPars);
            defaults.freqShift.sd      = 4*ones(1, nPars);
            defaults.freqShift.RegFun  = '';
                
            % Initialize metabolite amplitudes as free with value 0
            defaults.metAmpl.fun     = 'free';
            defaults.metAmpl.gradfun = 'free';
            defaults.metAmpl.lb      = 0.1*ones(1, nPars);
            defaults.metAmpl.ub      = 100*ones(1, nPars);
            defaults.metAmpl.init    = 0.25*ones(1, nPars);
            defaults.metAmpl.ex      = 0.5*ones(1, nPars);
            defaults.metAmpl.sd      = 100*ones(1, nPars);
            defaults.metAmpl.RegFun  = '';
    
            % Initialize baseline amplitudes as free with value 0
            defaults.baseAmpl.fun     = 'free';
            defaults.baseAmpl.gradfun = 'free';
            defaults.baseAmpl.lb      = -1*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.ub      = 1*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.init    = 0.1*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.ex      = 0.01*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.sd      = 100*ones(1, size(test{1}.BaselineBasis,2));
            defaults.baseAmpl.RegFun  = '';

            field_names = fieldnames(defaults.ph0);
            for ff = 1 : length(field_names)
                verifyEqual(testCase,test{1}.Options{1}.parametrizations.(parameter).(field_names{ff}), defaults.(parameter).(field_names{ff}) , ['Did set the default ' field_names{ff} ' from ' parameter ' correctly?'])
            end
        end

        function test_1D_external_parameter_per_basis_parsing(testCase,parameter_per_basis, fields)   
            
            % Find data and model procedure
            data = which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/data/1D/BSandMMs/NIfTIMRS/sim-001_ses-001_PRESS_3T_35_TE_A.nii.gz');
            model = which(['libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-model-procedures/ParameterParsing/1Step_invivo_external_' parameter_per_basis '.json']);
           
            parameter_name = strrep(parameter_per_basis,'_indiv','');
                     
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
            
            
            switch parameter_name
                case 'baseAmpl'
                    nPars = 11;
                    addition = linspace(0.01,0.01 * nPars,nPars);
                otherwise 
                    nPars = length(ModelProcedure.Steps{1}.basisset.include);
                    addition = linspace(0.01,0.01 * nPars,nPars);
            end
            
             % Set default initals
            % Initialize Lorentzian LB as constant with value
            defaults.lorentzLB.fun     = 'free';
            defaults.lorentzLB.gradfun = 'free';
            defaults.lorentzLB.lb      = ones(1, nPars);
            defaults.lorentzLB.ub      = 8*ones(1, nPars);
            defaults.lorentzLB.init    = 3.4*ones(1, nPars);
            defaults.lorentzLB.ex      = 3.4*ones(1, nPars);
            defaults.lorentzLB.sd      = 1000*ones(1, nPars);
            defaults.lorentzLB.RegFun  = '';
    
            % Initialize frequency shifts as constant with value 0 Hz
            defaults.freqShift.fun     = 'free';
            defaults.freqShift.gradfun = 'free';
            defaults.freqShift.lb      = -5*ones(1, nPars);
            defaults.freqShift.ub      = 5*ones(1, nPars);
            defaults.freqShift.init    = ones(1, nPars);
            defaults.freqShift.ex      = 2*ones(1, nPars);
            defaults.freqShift.sd      = 4*ones(1, nPars);
            defaults.freqShift.RegFun  = '';
                
            % Initialize metabolite amplitudes as free with value 0
            defaults.metAmpl.fun     = 'free';
            defaults.metAmpl.gradfun = 'free';
            defaults.metAmpl.lb      = zeros(1, nPars);
            defaults.metAmpl.ub      = 100*ones(1, nPars);
            defaults.metAmpl.init    = 0.25*ones(1, nPars);
            defaults.metAmpl.ex      = 0.5*ones(1, nPars);
            defaults.metAmpl.sd      = 100*ones(1, nPars);
            defaults.metAmpl.RegFun  = '';
    
            % Initialize baseline amplitudes as free with value 0
            defaults.baseAmpl.fun     = 'free';
            defaults.baseAmpl.gradfun = 'free';
            defaults.baseAmpl.lb      = -1*ones(1, nPars);
            defaults.baseAmpl.ub      = 1*ones(1, nPars);
            defaults.baseAmpl.init    = 0.1*ones(1, nPars);
            defaults.baseAmpl.ex      = 0.01*ones(1, nPars);
            defaults.baseAmpl.sd      = 100*ones(1, nPars);
            defaults.baseAmpl.RegFun  = '';


            % Update other fields
            if ~strcmp(fields, 'lb')
                ModelProcedure.Steps{1}.parametrizations.(parameter_name).lb = defaults.(parameter_name).lb(1);
                ModelProcedure.Steps{1}.parametrizations.(parameter_name).(fields) = ModelProcedure.Steps{1}.parametrizations.(parameter_name).(fields)*ones(1,nPars) + addition;
            end
            
            % Overwrite with expected individual values  
            if ~(strcmp(fields,'lb') && (strcmp(parameter_name,'freqShift') || strcmp(parameter_name,'baseAmpl')))
                defaults.(parameter_name).(fields) = defaults.(parameter_name).(fields) + addition;
            else
                defaults.(parameter_name).(fields) = defaults.(parameter_name).(fields) - addition;
            end
            
            % Run gLCM
            test = Osprey_gLCM({data},ModelProcedure,0,1);
        
            

            field_names = fieldnames(defaults.lorentzLB);
            for ff = 1 : length(field_names)
                verifyEqual(testCase,test{1}.Options{1}.parametrizations.(parameter_name).(field_names{ff}), defaults.(parameter_name).(field_names{ff}) , ['Did we parse the external value of ' field_names{ff} ' from ' parameter_name ' correctly?'],'AbsTol', 1e-10)
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