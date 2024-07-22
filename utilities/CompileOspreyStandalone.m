%% Compile Osprey
function CompileOspreyStandalone(OutputDir,SPM12Dir,CreateHBCD, CreateCMD, CreateGUI, PackageInstaller, WidgetsDir,GUILayoutDir)
%% CompileOspreyStandalone(OutputDir,SPM12Dir,CreateHBCD, CreateCMD, CreateGUI, PackageInstaller, WidgetsDir,GUILayoutDir)
%   This function compiles Osprey including a command window version, a
%   HBCD specific version with reduced basis sets, and a GUI version. The
%   compilation includes SPM12 for coregistration and segmentation. The GUI
%   version is only compiled if the path to the Widgets and the GUI Layout
%   Toolboxes is supplied. For the compilation to work you have to prepare
%   Osprey the following way:
% 
%   Remove the startup.m file from your Matlab folder as this interfers with the compilation and results in non-crashing errors in the compiled version.
%   Make sure to run a clear all before runnning the script.
% 
%   Add SPM12 to your path and include all subfolders!
% 
%   Update the version number.
% 
%
%   USAGE:
%       CompileOspreyStandalone(OutputDir,SPM12Dir,WidgetsDir,GUILayoutDir)
%
%   INPUTS:
%       OutputDir    = output directory.
%       SPM12Dir     = path to SPM12.
%       CreateHBCD   = Compile HBCD version
%       CreateCMD    = Compile CMD version
%       CreateGUI    = Compile GUI version   
%       PackageInstaller = Package installer for runtime download
%       WidgetsDir   = path to the Widgets Toolbox. Use the old version as
%       the newest version is not working with the compiler
%       GUILayoutDir = path to the GUI Layout Toolbox.
%
%   OUTPUTS:
%       Compiled Osprey.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2022-05-19)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2022-05-19: First version of the code.


%% 0. Setup folder strucutre
if nargin < 8
    GUILayoutDir = [];
    if nargin < 7
        WidgetsDir = [];
        if nargin < 6
            PackageInstaller = 1;
            if nargin < 5
                CreateGUI = 1;
                if nargin < 4
                    CreateCMD = 0;
                    if nargin < 3
                        CreateHBCD = 0;
                        if nargin < 2
                            error('ERROR: You need to indicate the output directory and SPM 12 path!!');
                        end
                    end
                end
            end
        end
    end
end


if ~exist(OutputDir)
    mkdir(OutputDir); 
end

if CreateHBCD
    OutputDirHBCD = fullfile(OutputDir,'HBCD');
    if exist(OutputDirHBCD,'Dir')
        rmdir(OutputDirHBCD,'s');
    end
    mkdir(OutputDirHBCD);
end

if CreateCMD
    OutputDirCmd = fullfile(OutputDir,'CommandLine');
    if exist(OutputDirCmd,'Dir')
        rmdir(OutputDirCmd,'s');
    end
    mkdir(OutputDirCmd);
end

if CreateGUI
    OutputDirGUI = fullfile(OutputDir,'GUI');
    if exist(OutputDirGUI,'Dir')
        rmdir(OutputDirGUI,'s');
    end
    mkdir(OutputDirGUI);
end

[SettingsDir,~,~] = fileparts(which('OspreySettings.m'));
AllDirs      = strsplit(SettingsDir, filesep);
OspreyDir       = strjoin(AllDirs(1:end-1), filesep); % parent folder (= Osprey directory)

%% 1. Prepare SPM12 for compilation
% The code is copied from spm_make_standalone.m
%-Input arguments
%--------------------------------------------------------------------------
contentsver = '';

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(SPM12Dir,'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
%-Get the list of toolbox directories
tbxdir = fullfile(SPM12Dir,'toolbox');
d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:numel(d)
    fi = spm_select('List',d{i},'.*_cfg_.*\.m$');
    if ~isempty(fi)
        ft = [ft(:); cellstr(fi)];
    end
end
%-Create code to insert toolbox config
if isempty(ft)
    ftstr = '';
else
    ft = spm_file(ft,'basename');
    ftstr = sprintf('%s ', ft{:});
end
fprintf(fid,'values = {%s};\n', ftstr);
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(SPM12Dir,'Contents.m'),...
               fullfile(SPM12Dir,'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end
if ~isempty(contentsver)
    % Format: 'xxxx (SPMx) dd-mmm-yyyy'
    f = fileread(fullfile(spm('Dir'),'Contents.txt'));
    f = regexprep(f,'% Version \S+ \S+ \S+',['% Version ' contentsver]);
    fid = fopen(fullfile(spm('Dir'),'Contents.txt'),'w');
    fprintf(fid,'%s',f);
    fclose(fid);
end

%==========================================================================
%-Trim FieldTrip
%==========================================================================
d = fullfile(SPM12Dir,'external','fieldtrip','compat');
d = cellstr(spm_select('FPList',d,'dir'));
for i=1:numel(d)
    f = spm_file(d{i},'basename');
    nrmv = strncmp(f,'matlablt',8);
    if nrmv
        [dummy,I] = sort({f(9:end),version('-release')});
        nrmv = I(1) == 2;
    end
    if ~nrmv
        [sts, msg] = rmdir(d{i},'s');
    end
end

% Copy the License files 
if CreateCMD
    copyfile(fullfile(OspreyDir,'LICENSE.md'),fullfile(OutputDirCmd,'OSPREY_LICENSE.md'));
    copyfile(fullfile(SPM12Dir,'LICENCE.txt'),fullfile(OutputDirCmd,'SPM12_LICENCE.txt'));
    copyfile(fullfile(SPM12Dir,'Contents.txt'),fullfile(OutputDirCmd,'SPM12_Contents.txt'));
end

if CreateHBCD
    copyfile(fullfile(OspreyDir,'LICENSE.md'),fullfile(OutputDirHBCD,'OSPREY_LICENSE.md'));
    copyfile(fullfile(SPM12Dir,'LICENCE.txt'),fullfile(OutputDirHBCD,'SPM12_LICENCE.txt'));
    copyfile(fullfile(SPM12Dir,'Contents.txt'),fullfile(OutputDirHBCD,'SPM12_Contents.txt'));
end

if CreateGUI
    copyfile(fullfile(OspreyDir,'LICENSE.md'),fullfile(OutputDirGUI,'OSPREY_LICENSE.md'));
    copyfile(fullfile(SPM12Dir,'LICENCE.txt'),fullfile(OutputDirGUI,'SPM12_LICENCE.txt'));
    copyfile(fullfile(SPM12Dir,'Contents.txt'),fullfile(OutputDirGUI,'SPM12_Contents.txt'));
end
%% 2. HBCD export
% Compile a HBCD specific version of Osprey (Copy basis set folder into application folder)
if CreateHBCD

appFile = fullfile(OspreyDir, 'utilities','RunOspreyJob.m');

opts = compiler.build.StandaloneApplicationOptions(appFile,...
    'OutputDir',OutputDirHBCD,...
    'ExecutableIcon',fullfile(OspreyDir, 'graphics','osprey.gif'),...
    'ExecutableSplashScreen',fullfile(OspreyDir, 'graphics','osprey.gif'),...
    'ExecutableVersion','2.6.2',...
    'ExecutableName','OspreyHBCD',...
    'AdditionalFiles',{ fullfile(SPM12Dir),...
                       fullfile(OspreyDir,'coreg'),...
                       fullfile(OspreyDir,'fit','code'),...
                       fullfile(OspreyDir,'graphics'),...
                       fullfile(OspreyDir,'job'),...
                       fullfile(OspreyDir,'libraries'),...
                       fullfile(OspreyDir,'load'),...
                       fullfile(OspreyDir,'overview'),...
                       fullfile(OspreyDir,'plot'),...
                       fullfile(OspreyDir,'process'),...
                       fullfile(OspreyDir,'quantify'),...
                       fullfile(OspreyDir,'seg'),...
                       fullfile(OspreyDir,'settings'),...
                       fullfile(OspreyDir,'utilities'),...
                       fullfile(OspreyDir,'CODE_OF_CONDUCT.md'),...
                       fullfile(OspreyDir,'LICENSE.md'),...
                       fullfile(OspreyDir,'README.md'),... 
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','compile_mex.m'),... We have to include the mex files manually
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.c'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexa64'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexmaci64'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexw64')},...  
    'AutoDetectDataFiles','On',...
    'TreatInputsAsNumeric','Off',...
    'Verbose','On');


    buildResults = compiler.build.standaloneApplication(opts);
end

%% 3. Command Line export
% Compile a command line version of Osprey
if CreateCMD
appFile = fullfile(OspreyDir, 'utilities','RunOspreyJob.m');

opts = compiler.build.StandaloneApplicationOptions(appFile,...
    'OutputDir',OutputDirCmd,...
    'ExecutableIcon',fullfile(OspreyDir, 'graphics','osprey.gif'),...
    'ExecutableSplashScreen',fullfile(OspreyDir, 'graphics','osprey.gif'),...
    'ExecutableVersion','2.6.2',...
    'ExecutableName','OspreyCMD',...
    'AdditionalFiles',{ fullfile(SPM12Dir),...
                       fullfile(OspreyDir,'coreg'),...
                       fullfile(OspreyDir,'fit','code'),...
                       fullfile(OspreyDir,'graphics'),...
                       fullfile(OspreyDir,'job'),...
                       fullfile(OspreyDir,'libraries'),...
                       fullfile(OspreyDir,'load'),...
                       fullfile(OspreyDir,'overview'),...
                       fullfile(OspreyDir,'plot'),...
                       fullfile(OspreyDir,'process'),...
                       fullfile(OspreyDir,'quantify'),...
                       fullfile(OspreyDir,'seg'),...
                       fullfile(OspreyDir,'settings'),...
                       fullfile(OspreyDir,'utilities'),...
                       fullfile(OspreyDir,'CODE_OF_CONDUCT.md'),...
                       fullfile(OspreyDir,'LICENSE.md'),...
                       fullfile(OspreyDir,'README.md'),... 
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','compile_mex.m'),... We have to include the mex files manually
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.c'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexa64'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexmaci64'),...
                       fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexw64')},...  
    'AutoDetectDataFiles','On',...
    'TreatInputsAsNumeric','Off',...
    'Verbose','On');


    buildResults = compiler.build.standaloneApplication (opts);

    if PackageInstaller
        compiler.package.installer(buildResults,'OutputDir',fullfile(OspreyDir,'OspreyCMD_installer'), 'InstallerName', 'OspreyCMD_install', 'RuntimeDelivery', 'installer',...
            'AuthorCompany','The Johns Hopkins University', 'Description', 'Installer for Osprey CMD version',...
            'InstallationNotes', 'Thank you for downloading Osprey. This installer will allow you to run the Osprey CMD as a standalone application. Please copy the "basissets" folder into the "application" folder after the installation.',...
            'Version', '2.6.2', 'ApplicationName','OspreyCMD')
    end
end



%% 4. GUI export
% Compile GUI Osprey which is only working if you supplied the path to the
% widgets and GUI Layout Toolboxes
if CreateGUI

    appFile = fullfile(OspreyDir, 'GUI','Osprey.m');
    opts = compiler.build.StandaloneApplicationOptions(appFile,...
        'OutputDir',OutputDirGUI,...
        'ExecutableIcon',fullfile(OspreyDir, 'graphics','osprey.gif'),...
        'ExecutableSplashScreen',fullfile(OspreyDir, 'graphics','osprey.gif'),...
        'ExecutableVersion','2.6.2',...
        'ExecutableName','Osprey',...
        'AdditionalFiles',{fullfile(WidgetsDir),...
                           fullfile(GUILayoutDir),...
                           fullfile(SPM12Dir),...
                           fullfile(OspreyDir,'GUI'),...
                           fullfile(OspreyDir,'coreg'),...
                           fullfile(OspreyDir,'fit','code'),...
                           fullfile(OspreyDir,'graphics'),...
                           fullfile(OspreyDir,'job'),...
                           fullfile(OspreyDir,'libraries'),...
                           fullfile(OspreyDir,'load'),...
                           fullfile(OspreyDir,'overview'),...
                           fullfile(OspreyDir,'plot'),...
                           fullfile(OspreyDir,'process'),...
                           fullfile(OspreyDir,'quantify'),...
                           fullfile(OspreyDir,'seg'),...
                           fullfile(OspreyDir,'settings'),...
                           fullfile(OspreyDir,'utilities'),...
                           fullfile(OspreyDir,'CODE_OF_CONDUCT.md'),...
                           fullfile(OspreyDir,'LICENSE.md'),...
                           fullfile(OspreyDir,'README.md'),... 
                           fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','compile_mex.m'),... We have to include the mex files manually
                           fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.c'),...
                           fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexa64'),...
                           fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexmaci64'),...
                           fullfile(OspreyDir,'libraries','L-BFGS-B-C','Matlab','lbfgsb_wrapper.mexw64')},...  
        'AutoDetectDataFiles','On',...
        'TreatInputsAsNumeric','Off',...
        'Verbose','On');

if isempty(GUILayoutDir) || isempty(WidgetsDir)
    error('ERROR: You need to indicate the oGUILayout and Widgets Toolbox directory to compile the GUI!!');
end


    buildResults = compiler.build.standaloneApplication(opts);

    if PackageInstaller
        compiler.package.installer(buildResults,'OutputDir',fullfile(OutputDirGUI,'OspreyGUI_installer'), 'InstallerName', 'OspreyGUI_install', 'RuntimeDelivery', 'web',...
            'AuthorCompany','The Johns Hopkins University', 'Description', 'Installer for Osprey GUI version',...
            'InstallationNotes', 'Thank you for downloading Osprey. This installer will allow you to run the Osprey GUI as a standalone application. Please copy the "basissets" folder into the "application" folder after the installation',...
            'Version', '2.6.2', 'ApplicationName','OspreyGUI')
    end
end


end