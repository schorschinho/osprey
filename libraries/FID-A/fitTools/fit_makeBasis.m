% fit_makeBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [BASIS] = fit_makeBasis(folder, addMMFlag, sequence, editTarget)
% 
% DESCRIPTION:
% Generates a basis set in FID-A structure. The code will search all *.mat
% files in the input folder for FID-A structures with simulated spectra. It
% also performs sanity checks on the simulation parameters, and returns
% warnings if parameters are not identical for all parameters.
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

function [BASIS] = fit_makeBasis(folder, addMMFlag, sequence, editTarget)

% Parse input arguments
if nargin < 4
    if strcmp(sequence, 'MEGA')
        editTarget = 'GABA';
    else
        editTarget = 'none';
    end
    if nargin < 3
        addMMFlag = 1;
        if nargin < 2
            sequence = 'unedited';
        end
    end
end

% Collect *.mat filenames from input folder
mat_files       = dir([folder filesep '*.mat']);
mat_filenames   = strcat(folder, filesep, {mat_files.name});
idx = contains(mat_filenames, 'Ref');
mat_filenames(idx) = [];
nMets           = length(mat_filenames);

% Loop over all *.mat filenames, load their data, store in a buffer
for kk = 1:nMets
    temp = load(mat_filenames{kk});
    % Multiplexed experiments (e.g. MEGA/HERMES) have more than one sub-basis 
    % simulated per metabolite. Find out how many:
    basisFct = fieldnames(temp);
    % Load the signals, DC-correct and store them in separate dimensions
    for ll = 1:length(basisFct)
        if isfield(temp.(basisFct{1}), 'centerFreq')
            buffer.centerFreq           = temp.(basisFct{1}).centerFreq;
        else
            temp.(basisFct{1}).centerFreq = 3;
            buffer.centerFreq = 3;
        end
        temp.(basisFct{ll}).specs      = fftshift(fft(temp.(basisFct{ll}).fids, [], 1), 1);
        spectralwidth = temp.(basisFct{ll}).spectralwidth;
        sz = temp.(basisFct{ll}).sz;
        Bo = temp.(basisFct{ll}).Bo;
        f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
        ppm=f/(Bo*42.577);
        ppm=ppm+4.68;
        temp.(basisFct{ll}).ppm = ppm - (4.68 - temp.(basisFct{1}).centerFreq);
        temp.(basisFct{ll})            = op_dccorr(temp.(basisFct{ll}),'p');
        buffer.fids(:,kk,ll)           = temp.(basisFct{ll}).fids;
        buffer.specs(:,kk,ll)          = temp.(basisFct{ll}).specs;
    end
    % The following field should always be the same for all sub-bases
    % (unless something is seriously flawed with the simulation code)
    buffer.t(:,kk)              = temp.(basisFct{1}).t;
    buffer.ppm(:,kk)            = temp.(basisFct{1}).ppm;
    buffer.spectralwidth(kk)    = temp.(basisFct{1}).spectralwidth;
    buffer.dwelltime(kk)        = temp.(basisFct{1}).dwelltime;
    buffer.n(kk)                = temp.(basisFct{1}).sz(1);
    buffer.linewidth(kk)        = temp.(basisFct{1}).linewidth;
    buffer.Bo(kk)               = temp.(basisFct{1}).Bo;
    buffer.seq{kk}              = temp.(basisFct{1}).seq;
    if isfield(temp.(basisFct{1}),'name')
        buffer.name{kk}             = temp.(basisFct{1}).name;
    else
        C = strsplit(mat_files(kk).name,'_');
        C = C{end};
        buffer.name{kk} = strrep(C,'.mat','');
    end
    buffer.te(kk)               = temp.(basisFct{1}).te;
    buffer.dims                 = temp.(basisFct{1}).dims;
    buffer.flags                = temp.(basisFct{1}).flags;
end

% Test whether parameters are the same across all basis functions; flag
% warning if they are not; write into basis set struct if they are.
seq_params = {'spectralwidth','dwelltime','n','linewidth','Bo','seq','te', 'centerFreq'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}));
    if length(unique_params) > 1
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params;
    end
end

% Test whether ppm and t aves are the same across all basis functions; flag
% error if they are not; write into basis set struct if they are.
seq_params = {'ppm','t'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}),'stable');
    if length(unique_params) ~= BASIS.n
        error('WARNING! One or more sequence parameters are not the same across all input basis functions.');
    else
        BASIS.(seq_params{pp}) = unique_params';
    end
end

% If chosen, add MM
if addMMFlag
    n = BASIS.n;
    sw = BASIS.spectralwidth;
    Bo = BASIS.Bo;
    centerFreq = BASIS.centerFreq;
    % The amplitude and FWHM values are determined as for the LCModel and
    % TARQUIN algorithms (see Wilson et al., MRM 2011).
    hzppm = Bo*42.577;
    
    % To scale the amplitudes correctly, we first need to determine the
    % area of the 3.027 ppm CH3 signal of creatine
    [CrArea] = detCrArea(buffer);
    oneProtonArea = CrArea/3;
    
    % Next, we determine the area of a Gaussian singlet with nominal area 1
    testGaussian    = op_gaussianPeak(n,sw,Bo,centerFreq,0.1*hzppm,0,1);
    testGaussian    = op_dccorr(testGaussian,'p');
    gaussianArea    = sum(real(testGaussian.specs));
    
    % Now we know the scaling factor to generate MM/lipid signals with the
    % correct relative scaling with respect to the CH3 signal
    MM09            = op_gaussianPeak(n,sw,Bo,centerFreq,0.14*hzppm,0.91,3*oneProtonArea/gaussianArea);
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
        buffer.name{nMets+rr}       = MMLips{rr};
        for qq = 1:length(basisFct)
            buffer.fids(:,nMets+rr,qq)   = MMBase.(MMLips{rr}).fids;
            buffer.specs(:,nMets+rr,qq)  = MMBase.(MMLips{rr}).specs;
        end
    end
    
    BASIS.flags.addedMM     = 1;
    BASIS.nMM               = length(MMLips);
    save_str = '_MM';
else
    BASIS.flags.addedMM     = 0;
    BASIS.nMM               = 0;
    save_str = '_noMM';
end

% If spectral editing has been performed, do the SUM and DIFF spectra here
if strcmp(sequence, 'MEGA')
    % Automatic recognition of on/off data based on NAA and water peaks
    switch editTarget
        case 'GABA'
            range_NAA = [1.9 2.1];
            pts_NAA = BASIS.ppm >= range_NAA(1) & BASIS.ppm <= range_NAA(end);
            idx_NAA = find(strcmp(buffer.name,'NAA'));
            maxA = max(real(buffer.specs(pts_NAA,idx_NAA,1)));
            maxB = max(real(buffer.specs(pts_NAA,idx_NAA,2)));
            if maxA/maxB < 0.1
                % this means A is on, B is off, indices need to be swapped
                switchOrder = [2 1];
            elseif maxA/maxB > 10
                % this means A is off, B is on, indices stay as the are
                switchOrder = [1 2];
            end
            % apply the switch order
            buffer.fids = buffer.fids(:,:,switchOrder);
            buffer.specs = buffer.specs(:,:,switchOrder);
        case 'GSH'
            range_w = [4.6 4.8];
            pts_w = BASIS.ppm >= range_w(1) & BASIS.ppm <= range_w(end);
            idx_w = find(strcmp(buffer.name,'H2O'));
            maxA = max(real(buffer.specs(pts_w,idx_w,1)));
            maxB = max(real(buffer.specs(pts_w,idx_w,2)));
            if maxA/maxB < 0.1
                % this means A is on, B is off, indices need to be swapped
                switchOrder = [2 1];
            elseif maxA/maxB > 10
                % this means A is off, B is on, indices stay as the are
                switchOrder = [1 2];
            end
            % apply the switch order
            buffer.fids = buffer.fids(:,:,switchOrder);
            buffer.specs = buffer.specs(:,:,switchOrder);
        case 'PE322'
            range_Cho = [3.0 3.4];
            pts_Cho = BASIS.ppm >= range_Cho(1) & BASIS.ppm <= range_Cho(end);
            idx_Cho = find(strcmp(buffer.name,'Cho'));
            maxA = max(real(buffer.specs(pts_Cho,idx_Cho,1)));
            maxB = max(real(buffer.specs(pts_Cho,idx_Cho,2)));
            if maxA/maxB < 0.1
                % this means A is on, B is off, indices need to be swapped
                switchOrder = [2 1];
            elseif maxA/maxB > 10
                % this means A is off, B is on, indices stay as the are
                switchOrder = [1 2];
            end
            % apply the switch order
            buffer.fids = buffer.fids(:,:,switchOrder);
            buffer.specs = buffer.specs(:,:,switchOrder);  
        case 'Lac'
            range_Lac = [4.0 4.3];
            pts_Lac = BASIS.ppm >= range_Lac(1) & BASIS.ppm <= range_Lac(end);
            idx_Lac = find(strcmp(buffer.name,'Lac'));
            maxA = max(real(buffer.specs(pts_Lac,idx_Lac,1)));
            maxB = max(real(buffer.specs(pts_Lac,idx_Lac,2)));
            if maxA/maxB < 0.1
                % this means A is on, B is off, indices need to be swapped
                switchOrder = [2 1];
            elseif maxA/maxB > 10
                % this means A is off, B is on, indices stay as the are
                switchOrder = [1 2];
            end
            % apply the switch order
            buffer.fids = buffer.fids(:,:,switchOrder);
            buffer.specs = buffer.specs(:,:,switchOrder);     
        case 'PE398'
            range_Ins = [3.8 4.1];
            pts_Ins = BASIS.ppm >= range_Ins(1) & BASIS.ppm <= range_Ins(end);
            idx_Ins = find(strcmp(buffer.name,'Ins'));
            maxA = max(real(buffer.specs(pts_Ins,idx_Ins,1)));
            maxB = max(real(buffer.specs(pts_Ins,idx_Ins,2)));
            if maxA/maxB < 0.1
                % this means A is on, B is off, indices need to be swapped
                switchOrder = [2 1];
            elseif maxA/maxB > 10
                % this means A is off, B is on, indices stay as the are
                switchOrder = [1 2];
            end
            % apply the switch order
            buffer.fids = buffer.fids(:,:,switchOrder);
            buffer.specs = buffer.specs(:,:,switchOrder);                      
    end
    
    % Now that we have guaranteed that the first dimension is always OFF
    % and the second one is always ON, we generate the DIFF and SUM.
    for rr = 1:length(buffer.name)
        switch editTarget
            case 'GABA'
                % GABA-edited data have co-edited MM09 that we need to put
                % in the DIFF. Therefore loop over metabolite names.
                % If MM09, then don't subtract the ON and OFF out; instead mimic
                % the co-edited MM09 signal in the DIFF by just a simple OFF MM09.
                MM09    = op_gaussianPeak(n,sw,Bo,centerFreq,0.085*hzppm,0.915,3*oneProtonArea/gaussianArea);
                MM09    = op_dccorr(MM09,'p');
                if strcmp(buffer.name{rr}, 'MM09')
                    buffer.fids(:,rr,3)      = MM09.fids; % DIFF
                    buffer.specs(:,rr,3)     = MM09.specs;
                else
                    buffer.fids(:,rr,3)      = buffer.fids(:,rr,2) - buffer.fids(:,rr,1); % DIFF
                    buffer.specs(:,rr,3)     = buffer.specs(:,rr,2) - buffer.specs(:,rr,1);
                end
            case 'GSH'
                % GSH-edited data have co-edited MM14 and MM12 that we need to put
                % in the DIFF. Therefore loop over metabolite names.
                % If one of those, then don't subtract the ON and OFF out; instead mimic
                % the co-edited signal in the DIFF by just a simple OFF MM14 or MM12.
                MM12     = op_gaussianPeak(n,sw,Bo,centerFreq,0.07*hzppm,1.20,2*oneProtonArea/gaussianArea);
                MM12     = op_dccorr(MM12,'p');
                MM14     = op_gaussianPeak(n,sw,Bo,centerFreq,0.095*hzppm,1.385,2*oneProtonArea/gaussianArea);
                MM14     = op_dccorr(MM14,'p');
                if strcmp(buffer.name{rr}, 'MM14')
                    buffer.fids(:,rr,3)      = MM14.fids; % DIFF
                    buffer.specs(:,rr,3)     = MM14.specs;
                elseif strcmp(buffer.name{rr}, 'MM12')
                    buffer.fids(:,rr,3)      = MM12.fids; % DIFF
                    buffer.specs(:,rr,3)     = MM12.specs;
                else
                    buffer.fids(:,rr,3)      = buffer.fids(:,rr,2) - buffer.fids(:,rr,1); % DIFF
                    buffer.specs(:,rr,3)     = buffer.specs(:,rr,2) - buffer.specs(:,rr,1);
                end
            case 'PE322'
                % PE-edited data 
                    buffer.fids(:,rr,3)      = buffer.fids(:,rr,2) - buffer.fids(:,rr,1); % DIFF
                    buffer.specs(:,rr,3)     = buffer.specs(:,rr,2) - buffer.specs(:,rr,1);   
            case 'PE398'
                % PE-edited data  have co-edited MM14 and MM12 that we need to put
                % in the DIFF. Therefore loop over metabolite names.
                % If one of those, then don't subtract the ON and OFF out; instead mimic
                % the co-edited signal in the DIFF by just a simple OFF MM14 or MM12.
                MM12     = op_gaussianPeak(n,sw,Bo,centerFreq,0.07*hzppm,1.20,2*oneProtonArea/gaussianArea);
                MM12     = op_dccorr(MM12,'p');
                MM14     = op_gaussianPeak(n,sw,Bo,centerFreq,0.095*hzppm,1.385,2*oneProtonArea/gaussianArea);
                MM14     = op_dccorr(MM14,'p');
                if strcmp(buffer.name{rr}, 'MM14')
                    buffer.fids(:,rr,3)      = MM14.fids; % DIFF
                    buffer.specs(:,rr,3)     = MM14.specs;
                elseif strcmp(buffer.name{rr}, 'MM12')
                    buffer.fids(:,rr,3)      = MM12.fids; % DIFF
                    buffer.specs(:,rr,3)     = MM12.specs;
                else
                    buffer.fids(:,rr,3)      = buffer.fids(:,rr,2) - buffer.fids(:,rr,1); % DIFF
                    buffer.specs(:,rr,3)     = buffer.specs(:,rr,2) - buffer.specs(:,rr,1);  
                end
            case 'Lac'
                % Lac-edited data dont have co-edited MMS that we need to put
                % in the DIFF. Therefore loop over metabolite names.
                    buffer.fids(:,rr,3)      = buffer.fids(:,rr,2) - buffer.fids(:,rr,1); % DIFF
                    buffer.specs(:,rr,3)     = buffer.specs(:,rr,2) - buffer.specs(:,rr,1);
        end
        buffer.fids(:,rr,4)      = buffer.fids(:,rr,2) + buffer.fids(:,rr,1); % SUM
        buffer.specs(:,rr,4)     = buffer.specs(:,rr,2) + buffer.specs(:,rr,1);
    end
elseif strcmp(sequence, 'HERMES') || strcmp(sequence, 'HERCULES')
    % Determine maximum signal intensities for water and NAA in each
    % sub-spectrum.
    range_NAA = [1.9 2.1];
    pts_NAA = BASIS.ppm >= range_NAA(1) & BASIS.ppm <= range_NAA(end);
    idx_NAA = find(strcmp(buffer.name,'NAA'));
    max_NAA(1) = max(real(buffer.specs(pts_NAA,idx_NAA,1)));
    max_NAA(2) = max(real(buffer.specs(pts_NAA,idx_NAA,2)));
    max_NAA(3) = max(real(buffer.specs(pts_NAA,idx_NAA,3)));
    max_NAA(4) = max(real(buffer.specs(pts_NAA,idx_NAA,4)));
    range_w = [4.6 4.8];
    pts_w = BASIS.ppm >= range_w(1) & BASIS.ppm <= range_w(end);
    idx_w = find(strcmp(buffer.name,'H2O'));
    max_w(1) = max(real(buffer.specs(pts_w,idx_w,1)));
    max_w(2) = max(real(buffer.specs(pts_w,idx_w,2)));
    max_w(3) = max(real(buffer.specs(pts_w,idx_w,3)));
    max_w(4) = max(real(buffer.specs(pts_w,idx_w,4)));
    % Sort the intensities in ascending order
    [~,order_w]   = sort(max_w);
    [~,order_NAA] = sort(max_NAA);
    % Now loop over the subspectra indices (A = 1, B = 2, etc) to determine
    % whether the respective experiments have high or low intensities:
    for ll = 1:4
        idx_w   = find(order_w == ll);
        idx_NAA = find(order_NAA == ll);
        
        if ismember(idx_w,[3 4])
            GSH.ON(ll) = 0;
        elseif ismember(idx_w,[1 2])
            GSH.ON(ll) = 1;
        end
        
        if ismember(idx_NAA,[3 4])
            GABA.ON(ll) = 0;
        elseif ismember(idx_NAA,[1 2])
            GABA.ON(ll) = 1;
        end
        
    end
    % Determine the sub-spectra indices belonging to each editing pattern
    idx_OFF_OFF = find(~GABA.ON & ~GSH.ON);
    idx_ON_OFF  = find(GABA.ON & ~GSH.ON);
    idx_OFF_ON  = find(~GABA.ON & GSH.ON);
    idx_ON_ON   = find(GABA.ON & GSH.ON);
    permVec = [idx_OFF_OFF, idx_ON_OFF, idx_OFF_ON, idx_ON_ON];

    % Commute the lines of fids and specs
    buffer.fids = buffer.fids(:,:,permVec);
    buffer.specs = buffer.specs(:,:,permVec);
    
    % Now generate the DIFFs and the SUM.
    % Making sure to include co-edited MMs appropriately in the DIFF.
    MM09    = op_gaussianPeak(n,sw,Bo,centerFreq,0.085*hzppm,0.915,3*oneProtonArea/gaussianArea);
    MM09    = op_dccorr(MM09,'p');
    MM12     = op_gaussianPeak(n,sw,Bo,centerFreq,0.07*hzppm,1.20,2*oneProtonArea/gaussianArea);
    MM12     = op_dccorr(MM12,'p');
    MM14     = op_gaussianPeak(n,sw,Bo,centerFreq,0.095*hzppm,1.385,2*oneProtonArea/gaussianArea);
    MM14     = op_dccorr(MM14,'p');
    for rr = 1:length(buffer.name)
        if strcmp(buffer.name{rr}, 'MM09')
            buffer.fids(:,rr,5)      = MM09.fids; % DIFF1 (GABA)
            buffer.specs(:,rr,5)     = MM09.specs;
            buffer.fids(:,rr,6)      = buffer.fids(:,rr,3) + buffer.fids(:,rr,4) - buffer.fids(:,rr,1) - buffer.fids(:,rr,2); % DIFF2 (GSH)
            buffer.specs(:,rr,6)     = buffer.specs(:,rr,3) + buffer.specs(:,rr,4) - buffer.specs(:,rr,1) - buffer.specs(:,rr,2);
        elseif strcmp(buffer.name{rr}, 'MM14')
            buffer.fids(:,rr,5)      = buffer.fids(:,rr,2) + buffer.fids(:,rr,4) - buffer.fids(:,rr,1) - buffer.fids(:,rr,3); % DIFF1 (GABA)
            buffer.specs(:,rr,5)     = buffer.specs(:,rr,2) + buffer.specs(:,rr,4) - buffer.specs(:,rr,1) - buffer.specs(:,rr,3);
            buffer.fids(:,rr,6)      = MM14.fids; % DIFF2 (GSH)
            buffer.specs(:,rr,6)     = MM14.specs;
        elseif strcmp(buffer.name{rr}, 'MM12')
            buffer.fids(:,rr,5)      = buffer.fids(:,rr,2) + buffer.fids(:,rr,4) - buffer.fids(:,rr,1) - buffer.fids(:,rr,3); % DIFF1 (GABA)
            buffer.specs(:,rr,5)     = buffer.specs(:,rr,2) + buffer.specs(:,rr,4) - buffer.specs(:,rr,1) - buffer.specs(:,rr,3);
            buffer.fids(:,rr,6)      = MM12.fids; % DIFF2 (GSH)
            buffer.specs(:,rr,6)     = MM12.specs;
        else
            buffer.fids(:,rr,5)      = buffer.fids(:,rr,2) + buffer.fids(:,rr,4) - buffer.fids(:,rr,1) - buffer.fids(:,rr,3); % DIFF1 (GABA)
            buffer.specs(:,rr,5)     = buffer.specs(:,rr,2) + buffer.specs(:,rr,4) - buffer.specs(:,rr,1) - buffer.specs(:,rr,3);
            buffer.fids(:,rr,6)      = buffer.fids(:,rr,3) + buffer.fids(:,rr,4) - buffer.fids(:,rr,1) - buffer.fids(:,rr,2); % DIFF2 (GSH)
            buffer.specs(:,rr,6)     = buffer.specs(:,rr,3) + buffer.specs(:,rr,4) - buffer.specs(:,rr,1) - buffer.specs(:,rr,2);
        end
        buffer.fids(:,rr,7)      = buffer.fids(:,rr,1) + buffer.fids(:,rr,3) + buffer.fids(:,rr,2) + buffer.fids(:,rr,4); % SUM
        buffer.specs(:,rr,7)     = buffer.specs(:,rr,1) + buffer.specs(:,rr,3) + buffer.specs(:,rr,2) + buffer.specs(:,rr,4);
    end
end

% Copy over the FID, specs, dims, and the metabolite names
BASIS.fids              = buffer.fids;
BASIS.specs             = buffer.specs;
BASIS.name              = buffer.name;
BASIS.dims              = buffer.dims;
BASIS.flags             = buffer.flags;
BASIS.nMets             = nMets;
BASIS.sz                = size(BASIS.fids);

% Normalize basis set
BASIS.scale = max(max(max(real(buffer.specs))));
BASIS.fids  = BASIS.fids ./ BASIS.scale;
BASIS.specs = BASIS.specs ./ BASIS.scale;

% Save as *.mat file
save(['BASIS' save_str '.mat'], 'BASIS');

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
ppm = in.ppm(:,1);
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