function[MRSCont] = osp_dqb_RunAll(MRSCont)
%% function[MRSCont] = osp_dqb_RunAll(MRSCont)
%
% Description: Function that applys the full data quality battery to all
% models called during Osprey Fitting. These include: 
%     1) "Vertical Shifts" on data and baseline
%     2) "Autocorrelation" on residual
%     3) "Runs test" on residual
%
% Input:     MRSCont = Osprey MRSCont with OspreyFit results
% Output:    MRSCont = Modified MRSCont with new field: "MRSCont.fit.dqb"
%
% Example usage:
%                   MRSCont = osp_dqb_RunAll(MRSCont);
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
MRSCont = [];
end

MRSCont.QM.dqb = []; % initialize empty struct

if MRSCont.flags.didFit
    MRSCont.QM.dqb.FitFlag = 1;
    switch MRSCont.opts.fit.method
        case 'Osprey_gLCM'
            for ms = 1:size(MRSCont.fit.results.metab,5)
                for ex = 1:size(MRSCont.fit.results.metab,4)
                    for ss = 1:size(MRSCont.fit.results.metab,3)
                        for kk = 1:size(MRSCont.fit.results.metab,2)
                            for mm = 1:size(MRSCont.fit.results.metab,1)
                                %% Ascertain signals of interest
                                Model = MRSCont.fit.results.metab{mm,kk,ss,ex,ms};
                                PPM = Model.Data.ppm;
                                Data = fftshift(fft(Model.Data.fids,[],1));
                                Residual = Model.Model{end}.fit.residual;
                                Baseline = Model.Model{end}.fit.baseline;
                                
                                % Select the relevant PPM range
                                RangeInd = [find(PPM>Model.Options{end}.optimFreqFitRange(1),1), find(PPM<Model.Options{end}.optimFreqFitRange(2),1,'last')];
                                %GapInd = ...;
                                if ~isempty(Model.Options{end}.gap)
                                    warning('Need to implement GAP in data quality metrics!')
                                end

                                %% LCM residual analytics
                                % Autocorrelation on residual
                                [MRSCont.QM.dqb.Res_SSAutocorr(mm,kk,ss,ex,ms),...
                                 MRSCont.QM.dqb.Res_MaxAutocorr(mm,kk,ss,ex,ms)] = osp_dqb_Autocorr(Residual, RangeInd);     
                                % "Runs" test on residual
                                MRSCont.QM.dqb.Residual_RunstestPVal(mm,kk,ss,ex,ms) = osp_dqb_RunsTest(Residual, RangeInd);
                                %% Nuisance signal analytics
                                % Calculates the residual-water-to-tCr ratio
                                MRSCont.QM.dqb.Water2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_SignalToResidualWater(Data, PPM);
                                % Calculates the Lipid-region-to-tCr ratio
                                MRSCont.QM.dqb.Lipid2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_SignalToLipid(Data, PPM);
                                %% LCM baseline analytics
                                % Calculates the baseline-to-tCr ratio
                                MRSCont.QM.dqb.Baseline2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_SignalToBaseline(Data, PPM, Baseline, RangeInd);
                                % Calculates the mean absolute curvature of the baseline
                                MRSCont.QM.dqb.MeanAbsCurvature(mm,kk,ss,ex,ms) = osp_dqb_BaselineCurvature(Baseline, PPM, RangeInd);
                                % Vertical shifts testing
                                [MRSCont.QM.dqb.anyNegative(mm,kk,ss,ex,ms),...
                                 MRSCont.QM.dqb.belowBaseline(mm,kk,ss,ex,ms)] = osp_dqb_Verticalshifts(Data,Baseline,RangeInd);
                            end
                        end
                    end
                end
            end
        otherwise
            % Only implemented the data quality battery for Osp_gLCM!!!
    end
else
    error('Must have modeled the data!')
end