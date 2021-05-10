function [MRSCont] = osp_fitWater(MRSCont, kk, fitWhich)
%% [MRSCont] = LCG_fitWater(MRSCont)
%   This function initializes and runs fitting of the water reference and
%   short-TE water scans.
%
%   USAGE:
%       [MRSCont] = osp_fitWater(MRSCont, fitWhich);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = Index of the cell containing the spectrum to be fit.
%       fitWhich    = String determining whether this function is performed
%                     on the water reference or the short-TE water scan.
%                     OPTIONS: - 'ref', 'w'
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-09)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-09: First version of the code.


% Find the basis set index corresponding to water
basisSet            = MRSCont.fit.basisSet;
h2o_idx = find(strcmp(basisSet.name, 'H2O'));
if isempty(h2o_idx)
    basisSet = load(MRSCont.opts.fit.basisSetFile);
    basisSet = basisSet.BASIS;
    h2o_idx = find(strcmp(basisSet.name, 'H2O'));
end
idx_toKeep = zeros(basisSet.nMets + basisSet.nMM,1);
idx_toKeep(h2o_idx) = 1;

% Remove the name, the FIDs and the specs of everything but water
% from the basis set.
basisSet.name   = basisSet.name(logical(idx_toKeep));
basisSet.fids   = basisSet.fids(:,logical(idx_toKeep),1); % index 1 because this is a GSH-OFF spectrum in edited data
basisSet.specs  = basisSet.specs(:,logical(idx_toKeep),1);
basisSet.nMets  = 1; basisSet.nMM = 0;

% Discern whether the reference or the short-TE water scan was selected
switch fitWhich
    case 'ref'
        str = 'ref';
    case 'w'
        str = 'w';
end

if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
%%  Construct the basis functions and the spectrum that is to be fit.
if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
else if MRSCont.flags.isMRSI == 1
        XVox = MRSCont.raw{1}.nXvoxels;
        YVox = MRSCont.raw{1}.nYvoxels;
        ZVox = MRSCont.raw{1}.nZvoxels;
    end
end

% Extract fit options
fitOpts     = MRSCont.opts.fit;
fitModel    = fitOpts.method;

if MRSCont.flags.isPRIAM == 1
    XVox = MRSCont.raw{1}.nXvoxels;
    for x = 1 : XVox
        % Apply scaling factor to the data
        dataToFit = op_takeVoxel(MRSCont.processed.(fitWhich){kk},x);      
        dataToFit = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
        % Call the fit function
        [fitParamsWater, resBasisSetWater]  = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);
        % Save back the fit parameters to MRSCont
        MRSCont.fit.resBasisSet{x}.(str).water{kk}      = resBasisSetWater;
        MRSCont.fit.results{x}.(str).fitParams{kk}   = fitParamsWater;
    end
elseif MRSCont.flags.isMRSI == 1
    fitOpts.isMRSI = 1;
    if isfield( MRSCont, 'mask')
        [r, c] = find(MRSCont.mask{1}== 1);
        cx = round(mean(r));
        cy = round(mean(c));
        cz = round(ZVox/2);
    else
       cx = round(XVox/2);
       cy = round(YVox/2);
       cz = round(ZVox/2);
    end
   if ZVox <=1
        dataToFit = op_takeVoxel(MRSCont.processed.(fitWhich){kk},[cx,cy]);  
    else
        dataToFit = op_takeVoxel(MRSCont.processed.(fitWhich){kk},[cx,cy,cz]); 
   end
    dataToFit = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
    [fitParamsWater, resBasisSetWater]  = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);
    if ZVox <=1
        MRSCont.fit.resBasisSet{cx,cy}.(str).water{kk}      = resBasisSetWater;
        MRSCont.fit.results{cx,cy}.(str).fitParams{kk}   = fitParamsWater;
    else  % 3D MRSI data
        MRSCont.fit.resBasisSet{cx,cy,cz}.(str).water{kk}      = resBasisSetWater;
        MRSCont.fit.results{cx,cy,cz}.(str).fitParams{kk}   = fitParamsWater;
    end
   
    % Fit all voxels
    for z = 1 : ZVox 
        for x = 1 : XVox 
            for y = 1 : YVox 
            [~] = printLog('OspreyFitWater',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
                try
                    for kk = 1 :MRSCont.nDatasets
                        if ZVox <=1
                            dataToFit = op_takeVoxel(MRSCont.processed.(fitWhich){kk},[x,y]);  
                        else
                            dataToFit = op_takeVoxel(MRSCont.processed.(fitWhich){kk},[x,y,z]); 
                        end
                    end
                    dataToFit = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                    [fitParamsWater, resBasisSetWater]  = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);

                    % Save back the fit parameters to MRSCont                
                    % 2D MRSI data
                    if ZVox <=1
                        MRSCont.fit.resBasisSet{x,y}.(str).water{kk}      = resBasisSetWater;
                        MRSCont.fit.results{x,y}.(str).fitParams{kk}   = fitParamsWater;
                    else  % 3D MRSI data
                        MRSCont.fit.resBasisSet{x,y,z}.(str).water{kk}      = resBasisSetWater;
                        MRSCont.fit.results{x,y,z}.(str).fitParams{kk}   = fitParamsWater;
                    end

                catch
                    if ZVox <=1
                        MRSCont.fit.results{x,y}.(str).fitParams{kk}   = MRSCont.fit.results{cx,cy}.(str).fitParams{kk};
                        MRSCont.fit.results{x,y}.(str).fitParams{kk}.ampl = 0;
                    else  % 3D MRSI data
                        MRSCont.fit.results{x,y,z}.(str).fitParams{kk}   = MRSCont.fit.results{cx,cy,cz}.(str).fitParams{kk};
                        MRSCont.fit.results{x,y,z}.(str).fitParams{kk}.ampl = 0;
                    end
                end
            end
        end
    end
 
else 
        % Apply scaling factor to the data
        dataToFit       = MRSCont.processed.(fitWhich){kk};
        dataToFit       = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});

        % Call the fit function
        [fitParamsWater, resBasisSetWater]  = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);
        % Write
        MRSCont.fit.resBasisSet.(str).water{kk}      = resBasisSetWater;
        MRSCont.fit.results.(str).fitParams{kk}   = fitParamsWater;
end
               
end