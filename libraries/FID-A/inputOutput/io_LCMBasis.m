% io_LCMBasis.m
% Helge Zollner, Johns Hopkins University 2020.
%
% USAGE:
% [BASIS] = io_LCMBasis(LCMBasisSet, addMMFlag, sequence, editTarget)
% 
% DESCRIPTION:
% Generates a basis set in FID-A structure. The code reads the LCM
% .BASIS-file and generates a Osprey compatible basis set. The output 
% BASIS file will be saved into the 'fit/basisset/user' folder. Next, you have
% to compare the field names of your BASIS file with the list in the fit_createMetabList()
% function ('libraries/FID-A/fitTools') and change them accordingly. Store 
% the updated  BASIS variable from your workspace. Afterwards,you have to change
% the path to your basisset file in the osp_fitInitialise() function in the '/fit' folder.
%
%
% INPUTS:
% folder    = folder containing *.mat files representing FID-A structures
% addMMFlag = Flag to decide whether MM and lipid basis functions should be
%               added to the basis set.
%             OPTIONS:  1 = Add MM+lip (Default)
%                       0 = Don't add MM+lip
% sequence  = sequence type
%             OPTIONS:  'unedited' (default)
%                       'MEGA'
%                       'HERMES'
%                       'HERCULES'
% editTarget= Target molecule of edited data.
%             OPTIONS:  'GABA'
%                       'GSH'
%                       '

%
% OUTPUTS:
% BASIS     = Simulated basis set in FID-A structure format. 

function [BASIS] = io_LCMBasis(LCMBasisSet, addMMFlag, sequence, editTarget, conjugate)
%Get Osprey folder
[settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
allFolders      = strsplit(settingsFolder, filesep);
ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)


% Parse input arguments
if nargin < 5
    conjugate = 1;
    if nargin < 4
        editTarget = 'none';
        if nargin < 3
            addMMFlag = 1;
            if nargin < 2
                sequence = 'unedited';
            end
        end
    end
end


     


% Load LCMBasisSet
[Read]=io_readlcmraw_basis(LCMBasisSet,conjugate);
disp('We are going to rename your basis set accroding to the Osprey naming convention now.\n')
Osp_names = ['Ala ','Asc ','Asp ', 'bHB ','bHG ', 'Cit ', 'Cr ', 'CrCH2 ', 'EtOH ', 'GABA \n',...
             'GPC ', 'GSH ', 'Glc ', 'Gln ', 'Glu ', 'Gly ', 'H2O ', 'mI ', 'Lac ', 'NAA ', 'NAAG \n', ...
             'PCh ', 'PCr ', 'PE ', 'Phenyl ', 'sI ', 'Ser ', 'Tau ', 'Tyros \n',...
             'MM09 ', 'MM12  ', 'MM14 ', 'MM17 ', 'MM20 ', 'Lip09 ', 'Lip13 ', 'Lip20 '];
filenames = fieldnames(Read);
for mm = 1 : length(fieldnames(Read))
     disp(['The name in your basis set is: ' filenames{mm}]);
     prompt = ['Type the name of the matching metabolite from the list or press return if no matchin metabolite is in the list: \n' Osp_names '\n'];
    name = input(prompt,'s');    
    if ~isempty(name)
        cRead.(name) = Read.(filenames{mm});
        ind = findstr(Osp_names,name);
        Osp_names(ind:(ind+length(name))) = [];
    end
end
filenames = fieldnames(cRead);
nMets     = length(fieldnames(cRead));

temp = cRead.(filenames{1});
BASIS.spectralwidth=temp.spectralwidth;
BASIS.dwelltime=1/temp.spectralwidth;
BASIS.n=temp.n;
BASIS.linewidth=temp.linewidth;
BASIS.Bo=temp.Bo;
BASIS.seq{1}=sequence;
BASIS.te=temp.te;
if isfield(temp, 'centerFreq')
    BASIS.centerFreq           = temp.centerFreq;
else
    BASIS.centerFreq           = 4.68;
end
BASIS.ppm=temp.ppm;
BASIS.t=temp.t;
BASIS.flags.writtentostruct=1;
BASIS.flags.gotparams=1;
BASIS.flags.leftshifted=1;
BASIS.flags.filtered=0;
BASIS.flags.zeropadded=0;
BASIS.flags.freqcorrected=0;
BASIS.flags.phasecorrected=0;
BASIS.flags.averaged=1;
BASIS.flags.addedrcvrs=1;
BASIS.flags.subtracted=1;
BASIS.flags.writtentotext=1;
BASIS.flags.downsampled=0;
BASIS.flags.isFourSteps=0;

    
% Loop over all *.mat filenames, load their data, store in a buffer
for kk = 1:nMets
    temp = cRead.(filenames{kk});
    temp.ppm = temp.ppm - (4.68 - BASIS.centerFreq);
    temp.centerFreq = BASIS.centerFreq;
    temp            = op_dccorr(temp,'p');
    BASIS.fids(:,kk)=temp.fids;
    BASIS.specs(:,kk)=temp.specs;
    BASIS.name{kk} = filenames{kk};

end
BASIS.ppm =temp.ppm;
BASIS.dims.t=1;
BASIS.dims.coils=0;
BASIS.dims.averages=0;
BASIS.dims.subSpecs=0;
BASIS.dims.extras=0;

BASIS.nMets = nMets;
BASIS.sz=[temp.sz(1) BASIS.nMets];

n = BASIS.n;
sw = BASIS.spectralwidth;
Bo = BASIS.Bo;
centerFreq = BASIS.centerFreq;
% The amplitude and FWHM values are determined as for the LCModel and
% TARQUIN algorithms (see Wilson et al., MRM 2011).
hzppm = Bo*42.577;

% To scale the amplitudes correctly, we first need to determine the
% area of the 3.027 ppm CH3 signal of creatine
[CrArea] = detCrArea(BASIS);
oneProtonArea = CrArea/3;

% Next, we determine the area of a Lorentzian singlet with nominal area 1
testLorentzian    = op_lorentzianPeak(n,sw,Bo,4.68,0.1*hzppm,0,1);
testLorentzian    = op_dccorr(testLorentzian,'p');
lorentzianArea    = sum(real(testLorentzian.specs));

% Find the creatine basis function
idx_H2O          = find(strcmp(BASIS.name,'H2O'));
if isempty(idx_H2O)
    H2O            = op_lorentzianPeak(n,sw,Bo,centerFreq,0.01*hzppm,4.68,2*oneProtonArea/lorentzianArea);
    H2O     = op_dccorr(H2O,'p');
    BASIS.name{nMets+1}       = 'H2O';
    BASIS.fids(:,nMets+1)   = H2O.fids;
    BASIS.specs(:,nMets+1)  = H2O.specs;
    nMets = nMets+1;
    BASIS.nMets = nMets;
    BASIS.sz=[temp.sz(1) BASIS.nMets];
end



% If chosen, add MM
if addMMFlag
    % Next, we determine the area of a Gaussian singlet with nominal area 1
    testGaussian    = op_gaussianPeak(n,sw,Bo,4.68,0.1*hzppm,0,1);
    testGaussian    = op_dccorr(testGaussian,'p');
    gaussianArea    = sum(real(testGaussian.specs));
    
    % Now we know the scaling factor to generate MM/lipid signals with the
    % correct relative scaling with respect to the CH3 signal
    MM09            = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm,0.91,3*oneProtonArea/gaussianArea);
%     MM09            = op_freqshift(MM09,-0.91*hzppm);
    MMBase.MM09     = op_dccorr(MM09,'p');
    MM12            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.21,2*oneProtonArea/gaussianArea);
    MMBase.MM12     = op_dccorr(MM12,'p');
    MM14            = op_gaussianPeak(n,sw,Bo,centerFreq,0.17*hzppm,1.43,2*oneProtonArea/gaussianArea);
    MMBase.MM14     = op_dccorr(MM14,'p');
    MM17            = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.67,2*oneProtonArea/gaussianArea);
    MMBase.MM17     = op_dccorr(MM17,'p');
    MM20a           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.08,1.33*oneProtonArea/gaussianArea);
    MM20b           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.25,0.33*oneProtonArea/gaussianArea);
    MM20c           = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.95,0.33*oneProtonArea/gaussianArea);
    MM20d           = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,3.0,0.4*oneProtonArea/gaussianArea);
    MM20            = op_addScans(MM20a,MM20b); MM20 = op_addScans(MM20,MM20c); MM20 = op_addScans(MM20,MM20d);
    MMBase.MM20     = op_dccorr(MM20,'p');
    Lip09           = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm,0.89,3*oneProtonArea/gaussianArea);
    MMBase.Lip09    = op_dccorr(Lip09,'p');
    Lip13a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,1.28,2*oneProtonArea/gaussianArea);
    Lip13b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.89*hzppm,1.28,2*oneProtonArea/gaussianArea);
    Lip13           = op_addScans(Lip13a,Lip13b);
    MMBase.Lip13    = op_dccorr(Lip13,'p');
    Lip20a          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.04,1.33*oneProtonArea/gaussianArea);
    Lip20b          = op_gaussianPeak(n,sw,Bo,centerFreq,0.15*hzppm,2.25,0.67*oneProtonArea/gaussianArea);
    Lip20c          = op_gaussianPeak(n,sw,Bo,centerFreq,0.2*hzppm,2.8,0.87*oneProtonArea/gaussianArea);
    Lip20           = op_addScans(Lip20a,Lip20b); Lip20 = op_addScans(Lip20,Lip20c);
    MMBase.Lip20    = op_dccorr(Lip20,'p');
    MMLips = {'MM09','MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20'};
    
    % Now copy over the names, fids, and specs into the basis set structure
    for rr = 1:length(MMLips)
        BASIS.name{nMets+rr}       = MMLips{rr};
        BASIS.fids(:,nMets+rr)   = MMBase.(MMLips{rr}).fids;
        BASIS.specs(:,nMets+rr)  = MMBase.(MMLips{rr}).specs;
    end
    
    BASIS.flags.addedMM     = 1;
    BASIS.nMM               = length(MMLips);
    save_str = '_MM';
else
    BASIS.flags.addedMM     = 0;
    BASIS.nMM               = 0;
    save_str = '_noMM';
end


BASIS.sz                = size(BASIS.fids);

% Normalize basis set
BASIS.scale = max(max(max(real(BASIS.specs))));
BASIS.fids  = BASIS.fids ./ BASIS.scale;
BASIS.specs = BASIS.specs ./ BASIS.scale;

%Reorder Fields
structorder = {'spectralwidth', 'dwelltime', 'n', 'linewidth', ...
                'Bo', 'seq', 'te', 'centerFreq', 'ppm', 't', 'flags', ...
                'nMM', 'fids', 'specs', 'name', 'dims','nMets', 'sz', 'scale'};

BASIS = orderfields(BASIS, structorder);

if ~exist(fullfile(ospFolder,'fit','basissets','user'),'dir')
    mkdir(fullfile(ospFolder,'fit','basissets','user'));
end
% Save as *.mat file
save(fullfile(ospFolder,'fit','basissets','user',['BASIS' save_str '.mat']), 'BASIS');

end


% detCrArea.m
% Georg Oeltzschner, Johns Hopkins University 2020
% 
% USAGE:
% [CrArea] = detCrArea(buffer);
% 
% DESCRIPTION:
% Finds the creatine spectrum in the temporary basis set buffer, then fits
% a Lorentzian to the 3.027 ppm CH3 creatine singlet to determine its area.
% Subsequently, macromolecule and lipid basis functions are scaled
% accordingly.
% 
% INPUTS:
% in        = a temporary buffer containing simulated basis functions
%
% OUTPUTS:
% CrArea    = Estimated area under the 3.027 ppm CH3 Cr singlet.


function [CrArea] = detCrArea(in);

% Find the creatine basis function
idx_Cr          = find(strcmp(in.name,'Cr'));
if isempty(idx_Cr)
    error('No basis function with nametag ''Cr'' found! Abort!');
end

%[~, idx_3027]   = min(abs(buffer.ppm(:,1)-3.027));

% Determine the window where we are going to look for the peak.
ppm = in.ppm(:);
ppmmin = 3.027 - 0.4;
ppmmax = 3.027 + 0.4;
refWindow = in.specs(ppm>ppmmin & ppm<ppmmax, idx_Cr);
ppmWindow = in.ppm(ppm>ppmmin & ppm<ppmmax);

% Find the maximum and its index
maxRef_index    = find(abs(real(refWindow)) == max(abs(real((refWindow)))));
maxRef          = real(refWindow(maxRef_index));

% Determine an initial estimate for the FWHM
% Peak lines can be super narrow, so overestimate it slightly
gtHalfMax   = find(abs(real(refWindow)) >= 0.4*abs(maxRef));
FWHM1       = abs(ppmWindow(gtHalfMax(1)) - ppmWindow(gtHalfMax(end)));
FWHM1       = FWHM1*(42.577*in.Bo(1));  %Assumes proton.

% Determine an initial estimate for the center frequency of the Cr peak
crFreq = ppmWindow(maxRef_index);

% Set up the fit
parsGuess=zeros(1,5);
parsGuess(1) = maxRef;  % amplitude
parsGuess(2) = (5*in.Bo/3)/(42.577*in.Bo); %FWHM.  Assumes Proton.  LW = 5/3 Hz/T.   % FWHM. Assumes Proton.
parsGuess(3) = crFreq;  % center frequency
parsGuess(4) = 0;       % baseline offset
parsGuess(5) = 0;       % phase
    
% Run first guess
yGuess  = op_lorentz(parsGuess, ppmWindow);
parsFit = nlinfit(ppmWindow, real(refWindow'), @op_lorentz, parsGuess);
yFit    = op_lorentz(parsFit, ppmWindow);
    
% figure;
% plot(ppmWindow,refWindow,'.',ppmWindow,yGuess,':',ppmWindow,yFit);
% legend('data','guess','fit');

CrArea = sum(yFit);

end