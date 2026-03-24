function MRSCont = osp_scale_yaxis(MRSCont,Module)
    switch Module
        case 'OspreyLoad' 
            ppmmax = 4.5;
            ppmmin = 0.2;   
            Range = zeros(2,MRSCont.nDatasets);
            for kk = 1 : MRSCont.nDatasets
                temp_spec = op_freqrange(MRSCont.raw{kk}, ppmmin, ppmmax);
                Range(1,kk) = op_findMaxMedian(temp_spec);
                Range(2,kk) = op_findMinMedian(temp_spec);           
            end
            MRSCont.plot.load.mets.max = Range(1,:);
            MRSCont.plot.load.mets.min = Range(2,:);
            MRSCont.plot.load.mets.ContMax = max(Range(1,:));
            MRSCont.plot.load.mets.ContMin = min(Range(2,:)); 

            if MRSCont.flags.hasRef
                ppmmax = 2*4.68;
                ppmmin = 0;   
                Range = zeros(2,MRSCont.nDatasets);
                for kk = 1 : MRSCont.nDatasets
                    temp_spec = op_freqrange(MRSCont.raw_ref{kk}, ppmmin, ppmmax);
                    Range(1,kk) = op_findMaxMedian(temp_spec);
                    Range(2,kk) = op_findMinMedian(temp_spec);            
                end
                MRSCont.plot.load.ref.max = Range(1,:);
                MRSCont.plot.load.ref.min = Range(2,:);
                MRSCont.plot.load.ref.ContMax = max(Range(1,:));
                MRSCont.plot.load.ref.ContMin = min(Range(2,:)); 
            end
            if MRSCont.flags.hasWater
                ppmmax = 2*4.68;
                ppmmin = 0;   
                Range = zeros(2,MRSCont.nDatasets);
                for kk = 1 : MRSCont.nDatasets
                    temp_spec = op_freqrange(MRSCont.raw_w{kk}, ppmmin, ppmmax);
                    Range(1,kk) = op_findMaxMedian(temp_spec,1);
                    Range(2,kk) = op_findMinMedian(temp_spec,1);             
                end
                MRSCont.plot.load.w.max = Range(1,:);
                MRSCont.plot.load.w.min = Range(2,:);
                MRSCont.plot.load.w.ContMax = max(Range(1,:));
                MRSCont.plot.load.w.ContMin = min(Range(2,:)); 
            end

        case 'OspreyProcess'
            SubSpecNames = fieldnames(MRSCont.processed);
            NoSubSpec = length(fieldnames(MRSCont.processed));
            for ss = 1 : NoSubSpec
                Range = zeros(2,MRSCont.nDatasets);
                switch SubSpecNames{ss}
                    case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum','mm'}
                        ppmmax = 4.5;
                    case {'ref', 'w'}
                        ppmmax = 2*4.68;        
                end
                switch SubSpecNames{ss}
                    case {'A', 'B', 'C', 'D', 'diff1', 'diff2', 'sum'}
                        ppmmin = 0.2;
                    case {'ref', 'w','mm'}
                        ppmmin = 0;        
                end
                for kk = 1 : MRSCont.nDatasets
                    temp_spec = op_freqrange(MRSCont.processed.(SubSpecNames{ss}){kk}, ppmmin, ppmmax);
                    if ~(strcmp(SubSpecNames{ss},'ref') || strcmp(SubSpecNames{ss},'w'))
                        Range(1,kk) = op_findMaxMedian(temp_spec);
                        Range(2,kk) = op_findMinMedian(temp_spec);
                    else
                        Range(1,kk) = op_findMaxMedian(temp_spec,1);
                        Range(2,kk) = op_findMinMedian(temp_spec,1);       
                    end
                end
                MRSCont.plot.processed.(SubSpecNames{ss}).max = Range(1,:);
                MRSCont.plot.processed.(SubSpecNames{ss}).min = Range(2,:);
                MRSCont.plot.processed.(SubSpecNames{ss}).ContMax = max(Range(1,:));
                MRSCont.plot.processed.(SubSpecNames{ss}).ContMin = min(Range(2,:)); 
            end
        case 'OspreyFit'

            MRSCont.plot.fit.match = 1; % Scaling between datasets is turned on by default
            if MRSCont.flags.isUnEdited
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'off');
            elseif MRSCont.flags.isMEGA
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'off');
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'diff1');
            elseif MRSCont.flags.isHERMES
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'sum');
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'diff1');
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'diff2');
            elseif MRSCont.flags.isHERCULES
                % For now, fit HERCULES like HERMES data
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'sum');
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'diff1');
                [MRSCont] = osp_extract_minmax_fit(MRSCont, 'diff2');
            end
            if MRSCont.flags.hasRef
               [MRSCont] = osp_extract_minmax_fit(MRSCont, 'ref'); 
            end
            if MRSCont.flags.hasWater
               [MRSCont] = osp_extract_minmax_fit(MRSCont, 'w'); 
            end
    end
end