function MRSCont = osp_scale_yaxis(MRSCont,Module)
    switch Module
        case 'OspreyLoad' 
            ppmmax = 4.5;
            ppmmin = 0.2;   
            Range = zeros(2,MRSCont.nDatasets(1),MRSCont.nDatasets(2));
            for kk = 1 : MRSCont.nDatasets(1)
                for ll = 1: 1:MRSCont.nDatasets(2)
                    metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
                    temp_spec = op_freqrange(MRSCont.raw{metab_ll,kk}, ppmmin, ppmmax);
                    Range(1,kk,metab_ll) = op_findMaxMedian(temp_spec);
                    Range(2,kk,metab_ll) = op_findMinMedian(temp_spec);      
                end
            end
            MRSCont.plot.load.mets.max = Range(1,:,:);
            MRSCont.plot.load.mets.min = Range(2,:,:);
            MRSCont.plot.load.mets.ContMax = max(max(Range(1,:)));
            MRSCont.plot.load.mets.ContMin = min(min(Range(2,:))); 

            if MRSCont.flags.hasRef
                ppmmax = 2*4.68;
                ppmmin = 0;   
                Range = zeros(2,MRSCont.nDatasets(1));
                for kk = 1 : MRSCont.nDatasets(1)
                    for ll = 1: 1:MRSCont.nDatasets(2)
                        ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                        temp_spec = op_freqrange(MRSCont.raw_ref{ref_ll,kk}, ppmmin, ppmmax);
                        Range(1,kk,ref_ll) = op_findMaxMedian(temp_spec);
                        Range(2,kk,ref_ll) = op_findMinMedian(temp_spec);    
                    end
                end
                MRSCont.plot.load.ref.max = Range(1,:,:);
                MRSCont.plot.load.ref.min = Range(2,:,:);
                MRSCont.plot.load.ref.ContMax = max(max(Range(1,:)));
                MRSCont.plot.load.ref.ContMin = min(min(Range(2,:))); 
            end
            if MRSCont.flags.hasWater
                ppmmax = 2*4.68;
                ppmmin = 0;   
                Range = zeros(2,MRSCont.nDatasets(1));
                for kk = 1 : MRSCont.nDatasets(1)
                    for ll = 1: 1:MRSCont.nDatasets(2)
                        w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                        temp_spec = op_freqrange(MRSCont.raw_w{w_ll,kk}, ppmmin, ppmmax);
                        Range(1,kk,w_ll) = op_findMaxMedian(temp_spec,1);
                        Range(2,kk,w_ll) = op_findMinMedian(temp_spec,1);     
                    end
                end
                MRSCont.plot.load.w.max = Range(1,:,:);
                MRSCont.plot.load.w.min = Range(2,:,:);
                MRSCont.plot.load.w.ContMax = max(max(Range(1,:)));
                MRSCont.plot.load.w.ContMin = min(min(Range(2,:))); 
            end
            if MRSCont.flags.hasMMRef
                ppmmax = 2*4.68;
                ppmmin = 0;   
                Range = zeros(2,MRSCont.nDatasets(1));
                for kk = 1 : MRSCont.nDatasets(1)
                    for ll = 1: 1:MRSCont.nDatasets(2)
                        mm_ref_ll = MRSCont.opts.MultipleSpectra.mm_ref(ll);
                        temp_spec = op_freqrange(MRSCont.raw_mm_ref{mm_ref_ll,kk}, ppmmin, ppmmax);
                        Range(1,kk,mm_ref_ll) = op_findMaxMedian(temp_spec,1);
                        Range(2,kk,mm_ref_ll) = op_findMinMedian(temp_spec,1);     
                    end
                end
                MRSCont.plot.load.mm_ref.max = Range(1,:,:);
                MRSCont.plot.load.mm_ref.min = Range(2,:,:);
                MRSCont.plot.load.mm_ref.ContMax = max(max(Range(1,:)));
                MRSCont.plot.load.mm_ref.ContMin = min(min(Range(2,:))); 
            end

        case 'OspreyProcess'
            SubSpecNames = fieldnames(MRSCont.processed);
            NoSubSpec = length(fieldnames(MRSCont.processed));
            for ss = 1 : NoSubSpec
                Range = zeros(2,MRSCont.nDatasets(1));
                switch SubSpecNames{ss}
                    case {'metab','mm'}
                        ppmmax = 4.5;
                    case {'ref', 'w', 'mm_ref'}
                        ppmmax = 2*4.68;        
                end
                switch SubSpecNames{ss}
                    case {'metab'}
                        ppmmin = 0.2;
                    case {'ref', 'w','mm','mm_ref'}
                        ppmmin = 0;        
                end
                for kk = 1 : MRSCont.nDatasets(1)
                    temp_spec = op_freqrange(MRSCont.processed.(SubSpecNames{ss}){kk}, ppmmin, ppmmax);
                    if ~(strcmp(SubSpecNames{ss},'ref') || strcmp(SubSpecNames{ss},'w')|| strcmp(SubSpecNames{ss},'mm_ref'))
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