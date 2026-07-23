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
                            if isempty(Model.Options{end}.gap)
                                GapInd = [];
                            else
                                GapInd = [find(PPM>Model.Options{end}.gap(1),1), find(PPM<Model.Options{end}.gap(2),1,'last')];
                            end
                            %GapInd = ...;
                            if ~isempty(Model.Options{end}.gap)
                                warning('Need to implement GAP in data quality metrics!')
                            end
                            
                            for dim = 1:size(Residual,2)
                                %% LCM residual analytics
                                % Autocorrelation on residual
                                [MRSCont.QM.dqb.Res_SSAutocorr(mm,kk,ss,ex,ms),...
                                 MRSCont.QM.dqb.Res_MaxAutocorr(mm,kk,ss,ex,ms)] = osp_dqb_Autocorr(Residual(:,dim), RangeInd, GapInd);     
                                % "Runs" test on residual
                                MRSCont.QM.dqb.Residual_RunstestPVal(mm,kk,ss,ex,ms) = osp_dqb_RunsTest(Residual(:,dim), RangeInd, GapInd);
                                %% Nuisance signal analytics
                                % Calculates the residual-water-to-tCr ratio
                                MRSCont.QM.dqb.Water2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_ResidualWater2tCr (Data(:,dim), PPM);
                                % Calculates the Lipid-region-to-tCr ratio
                                MRSCont.QM.dqb.Lipid2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_Lipid2tCr(Data(:,dim), PPM);
                                %% LCM baseline analytics
                                % Calculates the baseline-to-tCr ratio
                                MRSCont.QM.dqb.Baseline2tCr_Ratio(mm,kk,ss,ex,ms) = osp_dqb_Baseline2tCr(Data(:,dim), PPM, Baseline(:,dim), RangeInd);
                                % Calculates the mean absolute curvature of the baseline
                                MRSCont.QM.dqb.MeanAbsCurvature(mm,kk,ss,ex,ms) = osp_dqb_BaselineCurvature(Baseline(:,dim), PPM, RangeInd);
                                % Vertical shifts testing
                                [MRSCont.QM.dqb.anyNegative(mm,kk,ss,ex,ms),...
                                 MRSCont.QM.dqb.belowBaseline(mm,kk,ss,ex,ms)] = osp_dqb_Verticalshifts(Data(:,dim),Baseline(:,dim),RangeInd);
                            end
                        end
                    end
                end
            end
        end
    otherwise
        % Only implemented the data quality battery for Osp_gLCM!!!
end
end