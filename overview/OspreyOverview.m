function [MRSCont] = OspreyOverview(MRSCont)
%% [MRSCont] = OspreyOverview(MRSCont)
%   This function creates te data structre needed for the overview panel
%   in the GUI. It sorts the data by groups and performs the needed
%   statistics.
%
%   USAGE:
%       MRSCont = OspreyLoad(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2019-02-19)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-11-11: First version of the code.

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
    if nargin<1
        error('ERROR: no input Osprey container specified.  Aborting!!');
    end

%%% 2. INITIALIZE VARIABLES %%%
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
FitNames = fieldnames(MRSCont.fit.results);
NoFit = length(fieldnames(MRSCont.fit.results));

%%% 3. INTERPOLATION & NORMALIZATION %%%
% Processed data
MRSCont.overview.all_data = MRSCont.processed;
temp_sz = zeros(1,MRSCont.nDatasets);
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
        temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
    end
        [max_point,max_ind] = max(temp_sz);
    for kk = 1 : MRSCont.nDatasets
        if MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1) < max_point
            ppmRangeData        = MRSCont.processed.(SubSpecNames{ss}){1,max_ind}.ppm';
            ppmRangeDataToInt       = MRSCont.processed.(SubSpecNames{ss}){1,kk}.ppm;
            ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
            if sum(ppmIsInDataRange) == 0
                ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
            end
            MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
            MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.ppm = ppmRangeData;
            if mod(size(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),2)==0
                %disp('Length of vector is even.  Doing normal conversion');
                MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.fids=ifft(fftshift(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),[],MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t);
            else
                %disp('Length of vector is odd.  Doing circshift by 1');
                MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.fids=ifft(circshift(fftshift(MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.specs,MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t),1),[],MRSCont.overview.all_data.(SubSpecNames{ss}){1,kk}.dims.t);
            end

        end
    end
end
% Fits
for sf = 1 : NoFit
    for kk = 1 : MRSCont.nDatasets
        if strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w')
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}).water{kk};
        else if strcmp(FitNames{sf}, 'conc')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}){kk};
            else
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet.(FitNames{sf}){kk};
            end
        end
        fitParams   = MRSCont.fit.results.(FitNames{sf}).fitParams{kk};
        switch MRSCont.opts.fit.method
            case 'Osprey'
                ampl        = fitParams.ampl;
                zeroPhase   = fitParams.ph0;
                firstPhase  = fitParams.ph1;
                gaussLB     = fitParams.gaussLB;
                lorentzLB   = fitParams.lorentzLB;
                freqShift   = fitParams.freqShift;
                if strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w')
                    % If water, extract and apply nonlinear parameters
                    x       = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift];
                    [appliedBasisSet] = fit_applynlinwaterOsprey(basisSet, x, fitRangePPM);
                    fitted      = (appliedBasisSet.specs*ampl);
                    MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit      = fitted .* MRSCont.fit.scale{kk};
                    temp_spec     = op_freqrange(MRSCont.processed.(FitNames{sf}){1,kk}, fitRangePPM(1), fitRangePPM(end));
                    MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.ppm      = temp_spec.ppm;
                else if strcmp(FitNames{sf}, 'conc')
                        % If concatened metabolites, extract and apply nonlinear
                        % parameters and shifts
                        beta_j      = fitParams.beta_j;
                        spl_pos     = fitParams.spl_pos;
                        addFreqShift = fitParams.addFreqShift;
                        x          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j; addFreqShift];
                        [appliedBasisSet, B] = fit_applynlinOspreyMultiplex(basisSet, x, spl_pos, fitRangePPM);
                        % Convert to complex baseline and scaled basis functions
                        for mm = 1 : size(B,2)
                            B_conc = B(:,mm);
                            temp_appliedBasisSet.specs = appliedBasisSet.specs(:,:,mm);
                            temp_B      = reshape(B_conc, [length(B_conc)/2, 2]);
                            complex_B   = (temp_B(:,1) + 1i*temp_B(:,2));
                            fitted      = (temp_appliedBasisSet.specs*ampl + complex_B);
                            % Scale
                            complex_B   = complex_B .* MRSCont.fit.scale{kk};
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['fit' num2str(mm)])      = fitted .* MRSCont.fit.scale{kk};
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['baseline' num2str(mm)])      = real(complex_B);
                            temp_spec = op_freqrange(MRSCont.processed.A{1,kk}, fitRangePPM(1), fitRangePPM(end));
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['ppm' num2str(mm)])      = temp_spec.ppm;
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['res' num2str(mm)])      = real(temp_spec.specs) - MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['fit' num2str(mm)]);
                        end
                    else
                        % If metabolites, extract and apply nonlinear parameters
                        beta_j      = fitParams.beta_j;
                        spl_pos     = fitParams.spl_pos;
                        x          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j];
                        [appliedBasisSet, B] = fit_applynlinOsprey(basisSet, x, spl_pos, fitRangePPM);
                        % Convert to complex baseline and scaled basis functions
                        temp_B      = reshape(B, [length(B)/2, 2]);
                        complex_B   = (temp_B(:,1) + 1i*temp_B(:,2));
                        fitted      = (appliedBasisSet.specs*ampl + complex_B);
                        idx_1  = find(strcmp(appliedBasisSet.name,'NAA'));
                        idx_2  = find(strcmp(appliedBasisSet.name,'NAAG'));
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittNAA  = (ampl(idx_1).*appliedBasisSet.specs(:,idx_1).*MRSCont.fit.scale{kk})+(ampl(idx_2).*appliedBasisSet.specs(:,idx_2).*MRSCont.fit.scale{kk});
                        idx_1  = find(strcmp(appliedBasisSet.name,'Cr'));
                        idx_2  = find(strcmp(appliedBasisSet.name,'PCr'));
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittCr  = (ampl(idx_1).*appliedBasisSet.specs(:,idx_1).*MRSCont.fit.scale{kk})+(ampl(idx_2).*appliedBasisSet.specs(:,idx_2).*MRSCont.fit.scale{kk});                        
                        % Scale
                        complex_B   = complex_B .* MRSCont.fit.scale{kk};
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit      = fitted .* MRSCont.fit.scale{kk};
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.baseline      = real(complex_B);
                        temp_spec     = op_freqrange(MRSCont.processed.A{1,kk}, fitRangePPM(1), fitRangePPM(end));
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.ppm      = temp_spec.ppm;
                        specToPlot  = op_freqrange(MRSCont.processed.A{1,kk}, fitRangePPM(1), fitRangePPM(end));
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.res      = real(specToPlot.specs) - MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit;
                    end
                end
        end
    end
end

for sf = 1 : NoFit
    if strcmp(FitNames{sf}, 'conc')
        for kk = 1 : MRSCont.nDatasets
            temp_fit_sz.(FitNames{sf})(1,kk)= length(MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit1);
        end
        [max_point_fit.(FitNames{sf}),max_ind_fit.(FitNames{sf})] = max(temp_fit_sz.(FitNames{sf}));
    else
        for kk = 1 : MRSCont.nDatasets
            temp_fit_sz.(FitNames{sf})(1,kk)= length(MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit);
        end
        [max_point_fit.(FitNames{sf}),max_ind_fit.(FitNames{sf})] = max(temp_fit_sz.(FitNames{sf}));
    end
end

    
for sf = 1 : NoFit
    if strcmp(FitNames{sf}, 'conc')
            for kk = 1 : MRSCont.nDatasets
                if length(MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit1) < max_point_fit.(FitNames{sf})
                            ppmRangeData        = MRSCont.overview.all_fits.(FitNames{sf}){1,max_ind_fit.(FitNames{sf}).fit1}.ppm1';
                            ppmRangeDataToInt       = MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.ppm1;
                            ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                            if sum(ppmIsInDataRange) == 0
                                ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                            end
                            for mm = 1 : size(B,2)
                                MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['fit' num2str(mm)])      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['fit' num2str(mm)])(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                if ~strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w')
                                    MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['baseline' num2str(mm)]) = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['baseline' num2str(mm)])(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['res' num2str(mm)])      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['res' num2str(mm)])(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                                end
                                MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.(['ppm' num2str(mm)]) = ppmRangeData';
                            end
                end
            end
    else
        for kk = 1 : MRSCont.nDatasets
            if length(MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit) < max_point_fit.(FitNames{sf})
                        ppmRangeData        = MRSCont.overview.all_fits.(FitNames{sf}){1,max_ind_fit.(FitNames{sf})}.ppm';
                        ppmRangeDataToInt       = MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.ppm;
                        ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                        if sum(ppmIsInDataRange) == 0
                            ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                        end
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fit(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                        if ~strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w')
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.baseline = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.baseline(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittNAA      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittNAA(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittCr      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.fittCr(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                            MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.res      = interp1(ppmRangeDataToInt(ppmIsInDataRange), MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.res(ppmIsInDataRange), ppmRangeData, 'pchip', 'extrap');
                        end
                        MRSCont.overview.all_fits.(FitNames{sf}){1,kk}.ppm = ppmRangeData';
            end
        end
    end
end

MRSCont.overview.all_data_NAAnormalized = MRSCont.overview.all_data;
MRSCont.overview.all_fits_NAAnormalized = MRSCont.overview.all_fits;
for kk = 1 : MRSCont.nDatasets
    if isfield(MRSCont, 'quantify')
        if MRSCont.flags.isUnEdited
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.off{1,kk}.fittNAA= MRSCont.overview.all_fits.off{1,kk}.fittNAA/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.off{1,kk}.fittCr= MRSCont.overview.all_fits.off{1,kk}.fittCr/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.w{1,kk}.fit =  MRSCont.overview.all_fits.w{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA')); 
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                    MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                    MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                    MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                end
            end
            idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA')); 
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
            MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.fit= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
            MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.baseline= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
            MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.res= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
        end
        if MRSCont.flags.isMEGA
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.diff1{1,kk}.specs= MRSCont.overview.all_data.diff1{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.sum{1,kk}.specs= MRSCont.overview.all_data.sum{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.res1= MRSCont.overview.all_fits.conc{1,kk}.res1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.res2= MRSCont.overview.all_fits.conc{1,kk}.res2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                else
                    MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.res= MRSCont.overview.all_fits.diff1{1,kk}.res/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;                    
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;    
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.diff1{1,kk}.specs= MRSCont.overview.all_data.diff1{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.sum{1,kk}.specs= MRSCont.overview.all_data.sum{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.res1= MRSCont.overview.all_fits.conc{1,kk}.res1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.res2= MRSCont.overview.all_fits.conc{1,kk}.res2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                else
                    MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.res= MRSCont.overview.all_fits.diff1{1,kk}.res/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                       
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.diff1{1,kk}.specs= MRSCont.overview.all_data.diff1{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.sum{1,kk}.specs= MRSCont.overview.all_data.sum{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    if isfield(MRSCont.overview.all_fits, 'conc')
                        MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.res1= MRSCont.overview.all_fits.conc{1,kk}.res1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.res2= MRSCont.overview.all_fits.conc{1,kk}.res2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);                        
                    else
                        MRSCont.overview.all_fits.off{1,kk}.fit= MRSCont.overview.all_fits.off{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.off{1,kk}.baseline= MRSCont.overview.all_fits.off{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.off{1,kk}.res= MRSCont.overview.all_fits.off{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.diff1{1,kk}.res= MRSCont.overview.all_fits.diff1{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);                        
                    end
                end
            end
            if isfield(MRSCont.overview.all_fits, 'conc')
                idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
                MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_data_NAAnormalized.diff1{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.diff1{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_data_NAAnormalized.sum{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.sum{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.fit1= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.fit2= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.baseline1= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.baseline2= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.res1= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.res1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.res2= MRSCont.overview.all_fits_NAAnormalized.conc{1,kk}.res2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);                
            else
                idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
                MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_data_NAAnormalized.diff1{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.diff1{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_data_NAAnormalized.sum{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.sum{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.fit= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.fit= MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.baseline= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.baseline= MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.baseline/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.res= MRSCont.overview.all_fits_NAAnormalized.off{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.off(idx_1);
                MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.res= MRSCont.overview.all_fits_NAAnormalized.diff1{1,kk}.res/MRSCont.quantify.amplMets{1,kk}.off(idx_1);                
            end
        end
        if MRSCont.flags.isHERMES
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;                    
                else
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;                    
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;  
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                else
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl; 
                else
                    idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    if isfield(MRSCont.overview.all_fits, 'conc')
                        MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    else
                        MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                        MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                        MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                        MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                        MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                        MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                    end
                end
            end
            idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            if isfield(MRSCont.overview.all_fits, 'conc')
                MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            else
                MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
            end
        end
        if MRSCont.flags.isHERCULES
            if MRSCont.flags.hasRef
                MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                MRSCont.overview.all_data.ref{1,kk}.specs =  MRSCont.overview.all_data.ref{1,kk}.specs/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;                    
                else
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.baseline/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;                    
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.ref.fitParams{1,kk}.ampl;  
            else
                if MRSCont.flags.hasWater
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_data.w{1,kk}.specs =  MRSCont.overview.all_data.w{1,kk}.specs/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                if isfield(MRSCont.overview.all_fits, 'conc')
                    MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                else
                    MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                    MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.baseline/MRSCont.fit.results.w.fitParams{1,kk}.ampl;                    
                end
                MRSCont.overview.all_fits.ref{1,kk}.fit =  MRSCont.overview.all_fits.ref{1,kk}.fit/MRSCont.fit.results.w.fitParams{1,kk}.ampl;
                else
                    idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
                    MRSCont.overview.all_data.A{1,kk}.specs= MRSCont.overview.all_data.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.B{1,kk}.specs= MRSCont.overview.all_data.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.C{1,kk}.specs= MRSCont.overview.all_data.C{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    MRSCont.overview.all_data.D{1,kk}.specs= MRSCont.overview.all_data.D{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    if isfield(MRSCont.overview.all_fits, 'conc')
                        MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                        MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                    else
                        MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                        MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                        MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                        MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                        MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                        MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                    end
                end
            end
            idx_1  = find(strcmp(MRSCont.quantify.metabs,'tNAA'));
            MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.A{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.B{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.C{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs= MRSCont.overview.all_data_NAAnormalized.D{1,kk}.specs/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            if isfield(MRSCont.overview.all_fits, 'conc')
                MRSCont.overview.all_fits.conc{1,kk}.fit1= MRSCont.overview.all_fits.conc{1,kk}.fit1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.fit2= MRSCont.overview.all_fits.conc{1,kk}.fit2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.fit3= MRSCont.overview.all_fits.conc{1,kk}.fit3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline1= MRSCont.overview.all_fits.conc{1,kk}.baseline1/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline2= MRSCont.overview.all_fits.conc{1,kk}.baseline2/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
                MRSCont.overview.all_fits.conc{1,kk}.baseline3= MRSCont.overview.all_fits.conc{1,kk}.baseline3/MRSCont.quantify.amplMets{1,kk}.conc(idx_1);
            else
                MRSCont.overview.all_fits.diff1{1,kk}.fit= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                MRSCont.overview.all_fits.diff2{1,kk}.fit= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                MRSCont.overview.all_fits.sum{1,kk}.fit= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
                MRSCont.overview.all_fits.diff1{1,kk}.baseline= MRSCont.overview.all_fits.diff1{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff1(idx_1);
                MRSCont.overview.all_fits.diff2{1,kk}.baseline= MRSCont.overview.all_fits.diff2{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.diff2(idx_1);
                MRSCont.overview.all_fits.sum{1,kk}.baseline= MRSCont.overview.all_fits.sum{1,kk}.fit/MRSCont.quantify.amplMets{1,kk}.sum(idx_1);
            end
        end
    else
        error('This script works only on fully processed data. Run the whole Osprey pipeline first. Seg/Coreg is not needed')
    end
end

%%% 4. SORTING DATA  %%%
if MRSCont.flags.hasStatfile
    statCSV = readtable(MRSCont.file_stat, 'Delimiter', ',');
    group_idx = find(strcmp(statCSV{:,end},'group'));
    if isempty(group_idx)
        MRSCont.overview.groups = ones(MRSCont.nDatasets,1);
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    else
        MRSCont.overview.groups = statCSV{:,group_idx};
        MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
    end
else
    MRSCont.overview.groups = ones(MRSCont.nDatasets,1);
    MRSCont.overview.NoGroups = max(MRSCont.overview.groups);
end
MRSCont.overview.groupNames = cell(1,MRSCont.overview.NoGroups);
for g = 1 : MRSCont.overview.NoGroups
    MRSCont.overview.groupNames{g} = ['Group ' num2str(g)];
end

for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}) = MRSCont.overview.all_data.(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}) = MRSCont.overview.all_data_NAAnormalized.(SubSpecNames{ss})(1,MRSCont.overview.groups == g);
    end
end

for sf = 1 : NoFit
    for g = 1 : MRSCont.overview.NoGroups
        MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}) = MRSCont.overview.all_fits.(FitNames{sf})(1,MRSCont.overview.groups == g);
        MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}) = MRSCont.overview.all_fits_NAAnormalized.(FitNames{sf})(1,MRSCont.overview.groups == g);
    end
end

%%% 5. CALCULATE MEAN AND SD SPECTRA FOR VISUALIZATION %%%
%Data
for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        tempSubSpec = zeros(length(MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss})),MRSCont.overview.all_data.(SubSpecNames{1}){1,1}.sz(1));
        for kk = 1 : length(MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}))
          tempSubSpec(kk,:) = MRSCont.overview.(['sort_data_g' num2str(g)]).(SubSpecNames{ss}){1,kk}.specs;
        end
        MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' SubSpecNames{ss}]) = mean(real(tempSubSpec),1);
        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' SubSpecNames{ss}]) = std(real(tempSubSpec),1);
    end
MRSCont.overview.(['ppm_data_' SubSpecNames{ss}]) = MRSCont.overview.all_data.(SubSpecNames{ss}){1,1}.ppm;
end

for ss = 1 : NoSubSpec
    for g = 1 : MRSCont.overview.NoGroups
        tempSubSpec = zeros(length(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss})),MRSCont.overview.all_data_NAAnormalized.(SubSpecNames{1}){1,1}.sz(1));
        for kk = 1 : length(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}))
          tempSubSpec(kk,:) = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(SubSpecNames{ss}){1,kk}.specs;
        end
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' SubSpecNames{ss}]) = mean(real(tempSubSpec),1);
        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' SubSpecNames{ss}]) = std(real(tempSubSpec),1);
    end
end

%Fits
for sf = 1 : NoFit
    for g = 1 : MRSCont.overview.NoGroups
        if ~strcmp(FitNames{sf},'conc')
            tempSubSpec = zeros(length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1}),length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1}.ppm));
            tempSubBaseline = tempSubSpec;
            tempSubRes = tempSubSpec;
            for kk = 1 : length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}))
              tempSubSpec(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.fit;
              if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                tempSubBaseline(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.baseline;
                tempSubRes(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.res;
              end
            end
            MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' FitNames{sf}]) = mean(real(tempSubSpec),1);
            MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' FitNames{sf}]) = std(real(tempSubSpec),1);
            if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' FitNames{sf}]) = mean(real(tempSubBaseline),1);
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_baseline_' FitNames{sf}]) = std(real(tempSubBaseline),1);
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' FitNames{sf}]) = mean(real(tempSubRes),1);
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_res_' FitNames{sf}]) = std(real(tempSubRes),1);                 
            end
            MRSCont.overview.(['ppm_fit_' FitNames{sf}]) = MRSCont.overview.all_fits.(FitNames{sf}){1,1}.ppm;
        else
            for mm = 1 : length(fields(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1}))/3
                tempSubSpec = zeros(length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1}),length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1}.ppm1));
                tempSubBaseline = tempSubSpec;
                tempSubRes = tempSubSpec;
                for kk = 1 : length(MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}))
                  tempSubSpec(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.(['fit' num2str(mm)]);
                  if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                    tempSubBaseline(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.(['baseline' num2str(mm)]);
                    tempSubRes(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g)]).(FitNames{sf}){1,kk}.(['res' num2str(mm)]);
                  end
                end
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' FitNames{sf} '_fit' num2str(mm)]) = mean(real(tempSubSpec),1);
                MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' FitNames{sf} '_fit' num2str(mm)]) = std(real(tempSubSpec),1);
                if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                    MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' FitNames{sf} '_baseline_' num2str(mm)]) = mean(real(tempSubBaseline),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' FitNames{sf} '_baseline_' num2str(mm)]) = std(real(tempSubBaseline),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' FitNames{sf} '_res_' num2str(mm)]) = mean(real(tempSubRes),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g)]).(['sd_' FitNames{sf} '_res_' num2str(mm)]) = std(real(tempSubRes),1);
                end
            end
            MRSCont.overview.(['ppm_fit_' FitNames{sf}]) = MRSCont.overview.all_fits.(FitNames{sf}){1,1}.ppm1;
        end
    end
end

for sf = 1 : NoFit
    for g = 1 : MRSCont.overview.NoGroups
        if ~strcmp(FitNames{sf},'conc')
            tempSubSpec = zeros(length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1}),length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1}.ppm));
            tempSubBaseline = tempSubSpec;
            tempSubRes = tempSubSpec;
            for kk = 1 : length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}))
              tempSubSpec(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.fit;
              if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                tempSubBaseline(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.baseline;
                tempSubRes(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.res;
              end
            end
            MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_' FitNames{sf}]) = mean(real(tempSubSpec),1);
            MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_' FitNames{sf}]) = std(real(tempSubSpec),1);
            if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_baseline_' FitNames{sf}]) = mean(real(tempSubBaseline),1);
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_baseline_' FitNames{sf}]) = std(real(tempSubBaseline),1);
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_res_' FitNames{sf}]) = mean(real(tempSubRes),1);
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_res_' FitNames{sf}]) = std(real(tempSubRes),1);   
            end
            MRSCont.overview.(['ppm_fit_' FitNames{sf}]) = MRSCont.overview.all_fits.(FitNames{sf}){1,1}.ppm;
        else
            for mm = 1 : length(fields(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1}))/3
                tempSubSpec = zeros(length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1}),length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1}.ppm1));
                tempSubBaseline = tempSubSpec;
                tempSubRes = tempSubSpec;
                for kk = 1 : length(MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}))
                  tempSubSpec(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.(['fit' num2str(mm)]);
                  if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                    tempSubBaseline(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.(['baseline' num2str(mm)]);
                    tempSubRes(kk,:) = MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(FitNames{sf}){1,kk}.(['res' num2str(mm)]);
                  end
                end
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_' FitNames{sf} '_fit' num2str(mm)]) = mean(real(tempSubSpec),1);
                MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_' FitNames{sf} '_fit' num2str(mm)]) = std(real(tempSubSpec),1);
                if ~(strcmp(FitNames{sf}, 'ref') || strcmp(FitNames{sf}, 'w'))
                    MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_' FitNames{sf} '_baseline_' num2str(mm)]) = mean(real(tempSubBaseline),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_' FitNames{sf} '_baseline_' num2str(mm)]) = std(real(tempSubBaseline),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_' FitNames{sf} '_res_' num2str(mm)]) = mean(real(tempSubRes),1);
                    MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['sd_' FitNames{sf} '_res_' num2str(mm)]) = std(real(tempSubRes),1);                    
                end
            end
            MRSCont.overview.(['ppm_fit_' FitNames{sf}]) = MRSCont.overview.all_fits.(FitNames{sf}){1,1}.ppm1;
        end
    end

end
%%% 6. READ CORRELATION DATA INTO THE STRUCT %%%
if MRSCont.flags.hasStatfile
    cor = 1;
    while ~strcmp(statCSV{cor,end},'')
        MRSCont.overview.corr.Names{cor} = statCSV{cor,end};
        cor = cor + 1;
    end
    for cor = 1 : size(statCSV,2)-1
        MRSCont.overview.corr.Meas{cor} = statCSV{:,cor};
    end
end

%%% 7. CLEAN UP AND SAVE %%%
% Set exit flags
MRSCont.flags.didOverview          = 1;

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
save(fullfile(outputFolder, outputFile), 'MRSCont');

end
