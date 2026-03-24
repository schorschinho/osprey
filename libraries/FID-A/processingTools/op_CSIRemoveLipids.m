% op_CSIRemoveLipids.m
%
% Removes lipids from CSI data using L2 regularization. This minimizes the
% equation norm(x - x_0, 2) + beta * norm(W'x, 2)
%
% INPUT:
% MRSIStruct        = MRSI structure used in FID-A
% lipidComponents   = number of lipid spectra in the lipid basis
% lineWidthRange    = range of linewidth used for building lipid basis
% ppmRange          = ppm range used in building lipid basis
% beta              = regularization term
% plotBasis         = plot basis spectra
%                                                                                                                                                                                                                                                                          
%
% OUTPUT:
% MRSIStruct        = MRSI structure with lipids removed

function [MRSIStruct] = op_CSIRemoveLipids(MRSIStruct, basisArguments, lipidMask) 
   if nargin == 1
        basisArguments.Components = 2000;
        basisArguments.lineWidthRange = [10 80];
        basisArguments.PPMRange = [0.3 1.900];
        basisArguments.beta = 1;
        basisArguments.plotBasis = false;
        lipidMask = [];
   end
   if nargin == 2
        lipidMask = [];
   end
    % extract arguments from name value pairs
    lipidComponents = basisArguments.Components;
    lineWidthRange = basisArguments.lineWidthRange;
    lipidPPMRange = basisArguments.PPMRange;
    beta = basisArguments.beta;

    spectraSize = MRSIStruct.sz(MRSIStruct.dims.t);

    if isempty(lipidMask)
        % calculate lipid basis used for L2 regularization
        lipidBasis = createLipipBasis(MRSIStruct, lipidComponents, lineWidthRange, lipidPPMRange);
        if(basisArguments.plotBasis)
            figure
            plot(MRSIStruct.ppm, flip(real(lipidBasis),1));
        end
    
        % Not sure why I need this
        lipidBasis = cat(2,lipidBasis,flip(lipidBasis,1));
    else
        % Get lipid spectra from mask
        lipidBasisFids = reshape(MRSIStruct.fids,spectraSize,[]);
        lipidBasisSpecs = reshape(MRSIStruct.specs,spectraSize,[]);
        lipidMask  = squeeze(reshape(lipidMask,1,[]));
        lipidBasisFids = lipidBasisFids(:,lipidMask==1);
        lipidBasisSpecs = lipidBasisSpecs(:,lipidMask==1);
        lipidBasis = lipidBasisSpecs;
        temp = MRSIStruct;
        for j = 1 : size(lipidBasisFids,2)
            temp.fids = lipidBasisFids(:,j);
            temp.specs = lipidBasisSpecs(:,j);
            [~,spec_lip] = op_getLipid(temp,basisArguments.PPMRange,basisArguments.Components,floor(size(lipidBasisFids,1)*.75));
            lipidBasis(:,j)=spec_lip;
        end
        if(basisArguments.plotBasis)
            figure
            plot(MRSIStruct.ppm, abs(lipidBasis));
        end
    end

    % calculate solution from the basis
    L2Solution = inv(eye(spectraSize) + beta * (lipidBasis * lipidBasis'));

    data = reshape(MRSIStruct.specs,spectraSize,[]);
    data = L2Solution * data;
    if MRSIStruct.nZvoxels > 1
        data = reshape(data,spectraSize,MRSIStruct.sz(2),MRSIStruct.sz(3),MRSIStruct.sz(4));
    else if (MRSIStruct.nXvoxels * MRSIStruct.nYvoxels) > 1
            data = reshape(data,spectraSize,MRSIStruct.sz(2),MRSIStruct.sz(3));
        end
    end
    fids=ifft(fftshift(data,MRSIStruct.dims.t),[],MRSIStruct.dims.t);
    MRSIStruct.fids=fids;
    MRSIStruct.specs=data;

end

% calculate lipid basis
function lipidBasis = createLipipBasis(MRSIStruct, lipidComponents, lineWidthRange, lipidPPMRange)
    spectralWidth = MRSIStruct.spectralwidth;

    spectralPoints =  MRSIStruct.sz(MRSIStruct.dims.t);
    fidBasis = zeros(spectralPoints, lipidComponents);
    
    lipidStructure = load('Lip.mat', 'sysLip');
    lipidStructure = lipidStructure.sysLip;
    
    for iSpectra = 1:lipidComponents
        
        fidBasis(:, iSpectra) = getRandomLipidFids(spectralPoints, spectralWidth, ...
                                                  lineWidthRange, lipidPPMRange, ...
                                                  lipidStructure);
    end
    lipidBasis = fftshift(fft(fidBasis, [], 1), 1);
end

% calculate the single lipid spectra for the basis
function lipidFids = getRandomLipidFids(spectralPoints, spectralWidth, lineWidth, ...
                                        lipidPPMRange, lipidSystem)
    [randomLineWidth, randomPPM] = getRandomLineWidthandPPM(lineWidth, lipidPPMRange);
    lipidSystem.shifts = randomPPM;

    simulatedSignal = sim_onepulse(spectralPoints, spectralWidth, 3, randomLineWidth, lipidSystem);
    simulatedSignal = op_complexConj(simulatedSignal);
    simulatedSignal = addRandomPhase(simulatedSignal);
    simulatedSignal = scaleSpectra(simulatedSignal, randomPPM, lipidPPMRange);
    lipidFids = simulatedSignal.fids;
end


% pick a random number from lower bounds and upper bounds
function randomNumber = randomNumberInRange(lowerBounds, upperBounds)
    difference = upperBounds - lowerBounds;
    randomNumber = lowerBounds + difference * rand(1);
end
                
function simulatedSignal = addRandomPhase(simulatedSignal)
    sepctralPhase = randomNumberInRange(-180, 180);
    simulatedSignal = op_addphase(simulatedSignal, sepctralPhase, 0, 4.65, 1);
end

% scale spectra based on a normal distribution. Signal near the center of lipid
% range will be scaled high and signal near the edges scaled lower.
function simulatedSignal = scaleSpectra(simulatedSignal, lipidPPM, lipidRange)
    normalProbability = normpdf(lipidPPM, mean(lipidRange), diff(lipidRange)/4);
    simulatedSignal = op_ampScale(simulatedSignal, normalProbability);
    simulatedSignal = op_ampScale(simulatedSignal, 10);

end

function [randomLineWidth, randomPPM] = getRandomLineWidthandPPM(lineWidth, lipidPPMRange)
    randomLineWidth = randomNumberInRange(lineWidth(1), lineWidth(2));
    randomPPM = randomNumberInRange(lipidPPMRange(1), lipidPPMRange(2));
end
