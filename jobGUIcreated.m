%%% 0. CREDITS %%%
% This is a GUI generated Osprey jobFile. The code was kindly shared by Dr. Peter Van Schuerbeek (UZ Brussel)

%%% 1. SPECIFY SEQUENCE INFORMATION %%%
seqType = 'MEGA';
editTarget = {'GABA'};

%%% 2. SPECIFY DATA HANDLING AND MODELING OPTIONS %%%
opts.saveLCM = 0;
opts.savejMRUI = 0;
opts.saveVendor = 0;
opts.fit.method = 'Osprey';
opts.fit.includeMetabs = {'default'};
opts.fit.style = 'Separate';
opts.fit.range = [0.2 4.2];
opts.fit.rangeWater = [2.0 7.4];
opts.fit.bLineKnotSpace = 0.4;
opts.fit.fitMM = 1;

%%% 3. SPECIFY MRS DATA AND STRUCTURAL IMAGING FILES %%%
files = {'/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-01/mrs/sub-01_press/sub-01_PRESS30.dat',...
		 '/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-02/mrs/sub-02_press/sub-02_PRESS30.dat'};
files_ref = {'/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-01/mrs/sub-01_press/sub-01_PRESS30.dat',...
		 '/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-02/mrs/sub-02_press/sub-02_PRESS30.dat'};
files_w = {};
files_mm = {};
files_nii = {'/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-01/anat/sub-01_T1w.nii',...
		 '/Users/helge/Documents/GitHub/osprey/exampledata/twix/sub-02/anat/sub-02_T1w.nii'};

%%% 4. SPECIFY OUTPUT FOLDER %%%
outputFolder = '/Users/helge/Documents/GitHub/osprey/exampledata/twix/derivatives';