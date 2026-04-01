function [MRSCont] = OspreyMRSIHTMLReport(MRSCont,kk)
%% [MRSCont] = OspreyHTMLReport(MRSCont,kk)
%   This function creates a short HTML report of the processing and modeling
%   and should be called at the end of the analysis. It uses plotly to make
%   the results interactive.
%
%   USAGE:
%       MRSCont = OspreyHTMLReport(MRSCont,kk);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       kk          = subject index
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-10-31)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-10-31: First version of the code.
%% Prepartion
% Ensure that plotly works correctly
try
  if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
      VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                      round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                      round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
  else
      VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
  end

  out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
  PosLoadSpec = get(gcf,'Position');
  set(gcf,'Position',[PosLoadSpec(1) PosLoadSpec(2) 4*PosLoadSpec(4) PosLoadSpec(4)])
  p = fig2plotly(gcf, 'offline', true,'filename','Raw','fileopt','new','open',false);
catch
  close all;
  fprintf('Installing plotly for offline plotting.');
  getplotlyoffline('https://cdn.plot.ly/plotly-latest.min.js');
  try
    % Plot spectra
    if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
        VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                        round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                        round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
    else
        VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
    end

    out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
    PosLoadSpec = get(gcf,'Position');
    set(gcf,'Position',[PosLoadSpec(1) PosLoadSpec(2) 4*PosLoadSpec(4) PosLoadSpec(4)])
    p = fig2plotly(gcf, 'offline', true,'filename','Raw','fileopt','new','open',false);
    p = cleanup_spectra(p);
    movefile(fullfile(pwd,'Raw.html'),fullfile(outputFigures, 'Raw.html'));
    close(out)
  catch
    fprintf('Failed to install plotly for offline HTML plotting. Consult getplotlyoffline.');
  end
end

% Get colormaps and setup the inital path
colormaps = MRSCont.colormap;
ppmmin = 0.2;
ppmmax=4.2;

outputFolder    = fullfile(MRSCont.outputFolder,'Reports');
split_subject_path = strsplit(fileparts(MRSCont.files{kk}),filesep);
str_ind_sub = find(contains(split_subject_path,'sub'));
if ~isempty(str_ind_sub)
    sub_str = split_subject_path(str_ind_sub(1));
    sub_str = sub_str{1};
else
    sub_str = ['sub-' num2str(kk)];
end
outputFigures   = fullfile(MRSCont.outputFolder,'Reports','reportFigures',sub_str);
[foldername,filename,~]  = fileparts(MRSCont.files{kk});

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
if ~exist(outputFigures,'dir')
    mkdir(outputFigures);
end

% config = struct();
% config.responsive = true;
% config.fillFrame = true;
% config.frameMargins = 0;

%% OspreyLoad
if MRSCont.flags.didLoadData
    % Plot spectra
    if ~isfield(MRSCont.opts.MRSI.report,'VoxelIndices')
        VoxelIndices = [round(MRSCont.raw{kk}.nXvoxels/2),round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                        round(MRSCont.raw{kk}.nXvoxels/2)+1,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);
                        round(MRSCont.raw{kk}.nXvoxels/2)+2,round(MRSCont.raw{kk}.nYvoxels/2),round(MRSCont.raw{kk}.nZvoxels/2);];
    else
        VoxelIndices = MRSCont.opts.MRSI.report.VoxelIndices;
    end

    out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyLoad',1,'Fit1DStack',0,2,1,0,0);
    PosLoadSpec = get(gcf,'Position');
    set(gcf,'Position',[PosLoadSpec(1) PosLoadSpec(2) 4*PosLoadSpec(4) PosLoadSpec(4)])
    p = fig2plotly(gcf, 'offline', true,'filename','Raw','fileopt','new','open',false);
    p = cleanup_spectra(p);
    movefile(fullfile(pwd,'Raw.html'),fullfile(outputFigures, 'Raw.html'));
    close(out)

    % Quick Maps Load
    currentFolder = pwd;
    field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{1});
    for ff = 1 : length(field_names)
        MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{1}.(field_names{ff});
    end
    files_quickMaps_raw = [];
    names_quickMaps_raw = [];
    for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
        for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}))
            spec = MRSCont.opts.MRSI.quickMaps.specs{ss};
            name = MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll};
            names_quickMaps_raw{end+1} = [spec ' ' name];
            out = osp_plotQuickmaps(MRSCont, spec, name);
            p = fig2plotly(gcf, 'offline', true,'filename',[spec '_' name],'fileopt','new','open',false);
            p = cleanup_montages(p);
            movefile(fullfile(pwd,[spec '_' name '.html']),fullfile(outputFigures, [spec '_' name '.html']));
            files_quickMaps_raw{end+1} = fullfile(outputFigures, [spec '_' name '.html']);
            close(out)
        end
    end
end
%% OspreyCoreg
if MRSCont.flags.didCoreg
    out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 1, 0, 0, 0,1, MRSCont.raw{1, 1}.nZvoxels);
    exportgraphics(gcf, 'Coreg.png');
    CoregPos = get(gcf,'Position');
    % p = fig2plotly(gcf, 'offline', true,'filename','Coreg','fileopt','new','open',false);
    movefile(fullfile(pwd,'Coreg.png'),fullfile(outputFigures, 'Coreg.png'));
    close(out)
end
%% OspreySeg
if MRSCont.flags.didSeg
    out = osp_plotCoregMRSI(MRSCont, 'T1w_rMRSI', 0, 1, 1, 1,1, MRSCont.raw{1, 1}.nZvoxels);
    exportgraphics(gcf, 'CoregSeg.png');
    CoregSegPos = get(gcf,'Position');
    % p = fig2plotly(gcf, 'offline', true,'filename','CoregSeg','fileopt','new','open',false);
    movefile(fullfile(pwd,'CoregSeg.png'),fullfile(outputFigures, 'CoregSeg.png'));
    close(out)

    out = osp_plotSegMRSI(MRSCont,1,1,MRSCont.raw{1, 1}.nZvoxels, 1);
    set(gcf, 'Units', 'Normalized')
    SegPos = get(gcf,'OuterPosition');
    set(gcf, 'OuterPosition', [SegPos(1), SegPos(2), 1/5 + 0.05, 0.96])
    set(gcf, 'Units', 'Pixels')
    SegPos = get(gcf,'Position');
    p = fig2plotly(gcf, 'offline', true,'filename','Seg','fileopt','new','open',false);
    p = cleanup_seg(p);
    movefile(fullfile(pwd,'Seg.html'),fullfile(outputFigures, 'Seg.html'));
    close(out)
end
%% Osprey Process
if MRSCont.flags.didProcess
    % Plot spectra
    if isfield(MRSCont.processed,'A')
        TargetSpec = 'A';
    end
    if isfield(MRSCont.processed,'AFID')
        TargetSpec = 'AFID';
    end
    if isfield(MRSCont.processed,'diff1')
        TargetSpec = 'diff1';
    end
    out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyProcess',TargetSpec,'Fit1DStack',0,2,1,0,0);
    PosProcSpec = get(gcf,'Position');
    set(gcf,'Position',[PosProcSpec(1) PosProcSpec(2) 4*PosProcSpec(4) PosProcSpec(4)])
    p = fig2plotly(gcf, 'offline', true,'filename','Process','fileopt','new','open',false);
    p = cleanup_spectra(p);
    movefile(fullfile(pwd,'Process.html'),fullfile(outputFigures, 'Process.html'));
    close(out)

    % Quick Maps Process
    field_names = fieldnames(MRSCont.opts.MRSI.quickMapsList{2});
    for ff = 1 : length(field_names)
        MRSCont.opts.MRSI.quickMaps.(field_names{ff}) = MRSCont.opts.MRSI.quickMapsList{2}.(field_names{ff});
    end
    files_quickMaps_proc = [];
    names_quickMaps_proc = [];
    for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
        for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}))
            spec = MRSCont.opts.MRSI.quickMaps.specs{ss};
            name = MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll};
            names_quickMaps_proc{end+1} = [spec ' ' name];
            out = osp_plotQuickmaps(MRSCont, spec, name);
            p = fig2plotly(gcf, 'offline', true,'filename',[spec '_' name],'fileopt','new','open',false);
            p = cleanup_montages(p);
            movefile(fullfile(pwd,[spec '_' name '.html']),fullfile(outputFigures, [spec '_' name '.html']));
            files_quickMaps_proc{end+1} = fullfile(outputFigures, [spec '_' name '.html']);
            close(out)
        end
    end

    names_quickMaps_proc{end+1} = ['A SNR'];
    out = osp_plotQuickmaps(MRSCont, 'A', 'SNR');
    p = fig2plotly(gcf, 'offline', true,'filename',['A_SNR'],'fileopt','new','open',false);
    p = cleanup_montages(p);
    movefile(fullfile(pwd,'A_SNR.html'),fullfile(outputFigures,  'A_SNR.html'));
    files_quickMaps_proc{end+1} = fullfile(outputFigures, 'A_SNR.html');
    close(out)

    names_quickMaps_proc{end+1} = ['A FWHM'];
    out = osp_plotQuickmaps(MRSCont, 'A', 'FWHM');
    p = fig2plotly(gcf, 'offline', true,'filename',['A_FWHM'],'fileopt','new','open',false);
    p = cleanup_montages(p);
    movefile(fullfile(pwd,'A_FWHM.html'),fullfile(outputFigures,  'A_FWHM.html'));
    files_quickMaps_proc{end+1} = fullfile(outputFigures, 'A_FWHM.html');
    close(out)

    names_quickMaps_proc{end+1} = ['Global QC filtering'];
    out = osp_plotMetabolitemaps(MRSCont,'GlobalQC','tNAA_Acetyl_only',1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
    p = fig2plotly(gcf, 'offline', true,'filename',['Global_QC'],'fileopt','new','open',false);
    p = cleanup_montages(p);
    movefile(fullfile(pwd,'Global_QC.html'),fullfile(outputFigures,  'Global_QC.html'));
    files_quickMaps_proc{end+1} = fullfile(outputFigures, 'Global_QC.html');
    close(out)
end
%% OspreyFit
% Plot spectra
if MRSCont.flags.didFit
    out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices,'T1w_rMRSI','OspreyFit','metab','Fit1DStack',0,2,1,0,0);
    PosFitSpec = get(gcf,'Position');
    set(gcf,'Position',[PosFitSpec(1) PosFitSpec(2) 4*PosFitSpec(4) PosFitSpec(4)])
    p = fig2plotly(gcf, 'offline', true,'filename','Fit','fileopt','new','open',false);
    p = cleanup_spectra(p);
    movefile(fullfile(pwd,'Fit.html'),fullfile(outputFigures, 'Fit.html'));
    close(out)
end
%% OspreyQuantify
if MRSCont.flags.didQuantify
    if ~isfield(MRSCont.opts.MRSI.report, 'quantifications')
        quantifcations = {'tCr'};
    else
        quantifcations = MRSCont.opts.MRSI.report.quantifications;
    end

    if ~isfield(MRSCont.opts.MRSI.report, 'metabolites')
        metabolites = {'tNAA'};
    else
        metabolites = MRSCont.opts.MRSI.report.metabolites;
    end

    files_Quantification_Maps=[];
    names_Quantification_Maps=[];
    % Loop over quantifications
    for qq = 1 : length(quantifcations)
        for mm = 1 : length(metabolites)
            out = osp_plotMetabolitemaps(MRSCont,quantifcations{qq},metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
            p = fig2plotly(gcf, 'offline', true,'filename',[quantifcations{qq},'_',metabolites{mm}],'fileopt','new','open',false);
            p = cleanup_montages(p);
            movefile(fullfile(pwd,[quantifcations{qq},'_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_',metabolites{mm} '.html']));
            files_Quantification_Maps{end+1} = fullfile(outputFigures, [quantifcations{qq},'_',metabolites{mm} '.html']);
            names_Quantification_Maps{end+1} = [quantifcations{qq},' ',metabolites{mm}];
            close(out)
        end
    end

    files_Quantification_Maps_QC=[];
    names_Quantification_Maps_QC=[];
    % Loop over quantifications
    for qq = 1 : length(quantifcations)
        for mm = 1 : length(metabolites)
            out = osp_plotMetabolitemaps(MRSCont,[quantifcations{qq} '_QC'],metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
            p = fig2plotly(gcf, 'offline', true,'filename',[quantifcations{qq},'_QC','_',metabolites{mm}],'fileopt','new','open',false);
            p = cleanup_montages(p);
            movefile(fullfile(pwd,[quantifcations{qq},'_QC','_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_QC','_',metabolites{mm} '.html']));
            files_Quantification_Maps{end+1} = fullfile(outputFigures, [quantifcations{qq},'_QC','_',metabolites{mm} '.html']);
            names_Quantification_Maps{end+1} = [quantifcations{qq},' ',metabolites{mm}, 'QC filter'];
            close(out)
        end
    end

    files_Quantification_Maps_QCfilt=[];
    names_Quantification_Maps_QCfilt=[];
    % Loop over quantifications
    for qq = 1 : length(quantifcations)
        for mm = 1 : length(metabolites)
            out = osp_plotMetabolitemaps(MRSCont,[quantifcations{qq} '_QCfilt'],metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
            p = fig2plotly(gcf, 'offline', true,'filename',[quantifcations{qq},'_QCfilt','_',metabolites{mm}],'fileopt','new','open',false);
            p = cleanup_montages(p);
            movefile(fullfile(pwd,[quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']),fullfile(outputFigures,  [quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']));
            files_Quantification_Maps_QCfilt{end+1} = fullfile(outputFigures, [quantifcations{qq},'_QCfilt','_',metabolites{mm} '.html']);
            names_Quantification_Maps_QCfilt{end+1} = [quantifcations{qq},' ',metabolites{mm}, 'QC filtered'];
            close(out)
        end
    end

    % CRLB maps
    files_CRLB_Maps=[];
    names_CRLB_Maps=[];
    for mm = 1 : length(metabolites)
        out = osp_plotMetabolitemaps(MRSCont,'CRLBs',metabolites{mm},1,MRSCont.raw{1, 1}.nZvoxels,1,1,0);
        p = fig2plotly(gcf, 'offline', true,'filename',['CRLBs','_',metabolites{mm}],'fileopt','new','open',false);
        p = cleanup_montages(p);
        movefile(fullfile(pwd,['CRLBs','_',metabolites{mm} '.html']),fullfile(outputFigures,  ['CRLBs','_',metabolites{mm} '.html']));
        files_CRLB_Maps{end+1} = fullfile(outputFigures, ['CRLBs','_',metabolites{mm} '.html']);
        names_CRLB_Maps{end+1} = ['CRLBs', ' ',  quantifcations{qq},' ',metabolites{mm}, 'QC filtered'];
        close(out)
    end
end
%% OspreyOverview
if MRSCont.flags.didOverview

    % Global Concentration Results
    if isfield(MRSCont.opts.MRSI,'GlobalConc')
        quantifcations = MRSCont.opts.MRSI.GlobalConc.quantities;
        metabolites = MRSCont.opts.MRSI.GlobalConc.metabolites;
        files_GlobalConc=[];
        names_GlobalConc=[];
        for qq = 1 : length(quantifcations)
            for mm = 1 : length(metabolites)
                out = osp_plotGlobalConcentration(MRSCont,metabolites{mm},quantifcations{qq});
                set(gcf, 'Units', 'Normalized')
                Pos = get(gcf,'OuterPosition');
                set(gcf, 'OuterPosition', [Pos(1), Pos(2),0.4, 0.4])
                set(gcf, 'Units', 'Pixels')
                PosGlobalConc = get(gcf,'Position');
                exportgraphics(gcf, ['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']);
                 movefile(fullfile(pwd,['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']),fullfile(outputFigures,['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']));
                files_GlobalConc{end+1} = fullfile(outputFigures, ['GlobalConc_' quantifcations{qq},'_',metabolites{mm} '.png']);
                names_GlobalConc{end+1} = ['Global Concentration ', quantifcations{qq},' ',metabolites{mm}];
                close(out)
            end
        end
    end

    % Atlas analysis results
    if isfield(MRSCont.opts.MRSI,'atlas')
        quantifcations = MRSCont.opts.MRSI.atlas.quantities;
        regions = MRSCont.opts.MRSI.report.atlasregion;
        files_atlas=[];
        names_atlas=[];
        for qq = 1 : length(quantifcations)
            for rr = 1 : length(regions)
                out = osp_plotInteractiveAtlasAnalysis(MRSCont,1,quantifcations{qq},MRSCont.opts.MRSI.report.metabolites,'T1w_rMRSI',regions{rr});
                set(gcf, 'Units', 'Normalized')
                Pos = get(gcf,'OuterPosition');
                set(gcf, 'OuterPosition', [Pos(1), Pos(2),0.4, 0.4])
                set(gcf, 'Units', 'Pixels')
                PosAtlas = get(gcf,'Position');
                 exportgraphics(gcf, ['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']);
                movefile(fullfile(pwd,['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']),fullfile(outputFigures,['AtlasResults_' quantifcations{qq},'_region_',num2str(rr),'.png']));
                files_atlas{end+1} = fullfile(outputFigures, ['AtlasResults_' quantifcations{qq},'_region_',num2str(rr), '.png']);
                names_atlas{end+1} = ['Atlas Results ', quantifcations{qq},' ',regions{rr}];
                close(out)
            end
        end
    end
end
%% Write report in HTML files
%Write as relative path
outputFigures   = fullfile('reportFigures',sub_str);
%write an html report:
fid=fopen(fullfile(outputFolder,[sub_str,'-report.html']),'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
fprintf(fid,'\n<body>');

logoPath=which('osprey.png');
if ~isempty(logoPath)
    fprintf(fid,'\n<img src= " %s " width="35" height="28"> <b> \tOsprey MRSI Analysis Report</b> ',logoPath);
else
    fprintf(fid,'\n<b> Osprey MRSI Analysis Report</b>');
end
fprintf(fid,'\n<p><b>DATE:</b> %s \t <b>FILENAME:</b> %s </p>',date,filename);

if MRSCont.flags.didLoadData
    % OspreyLoad
    fprintf(fid,'\n<h2> Osprey Load</h2>');
    fprintf(fid,'\n<h3> Example Raw Spectra </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'Raw.html'),PosLoadSpec(4)*1.5);

    fprintf(fid,'\n<h3> Quickmaps Raw Data (Amplitude Integration)</h3>');
    for ff = 1 : length(files_quickMaps_raw)
        fprintf(fid,'\n %s ', names_quickMaps_raw{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_quickMaps_raw{ff});
    end
end

if MRSCont.flags.didCoreg
    % OspreyCoreg
    fprintf(fid,'\n<h2> Osprey Coregistration</h2>');
    fprintf(fid,'\n<h3> MRSI slice localization </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'Coreg.png'),CoregPos(4)*1.8);
end

if MRSCont.flags.didSeg
    % OspreySeg
    fprintf(fid,'\n<h2> Osprey Segmentation</h2>');
    fprintf(fid,'\n<h3> Outer mask + automated brain mask </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'CoregSeg.png'),CoregSegPos(4)*1.8);

    fprintf(fid,'\n<h3> Tissue fraction maps + automated masks </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'Seg.html'),SegPos(4)*1.5);
end

if MRSCont.flags.didProcess
    % OspreyProcess
    fprintf(fid,'\n<h2> Osprey Process</h2>');
    fprintf(fid,'\n<h3> Example Processed Spectra %s</h3>',TargetSpec);
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'Process.html'),PosProcSpec(4)*1.5);

    fprintf(fid,'\n<h3> Quickmaps Processed Data (Amplitude Integration)</h3>');
    for ff = 1 : length(files_quickMaps_proc)
        if ff ==  length(files_quickMaps_proc)-2
            fprintf(fid,'\n<h3> Quality Metric Maps </h3>');
        end
        fprintf(fid,'\n %s ', names_quickMaps_proc{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_quickMaps_proc{ff});
    end
end

if MRSCont.flags.didFit
    % OspreyFit
    fprintf(fid,'\n<h2> Osprey Fit</h2>');
    fprintf(fid,'\n<h3> Example Fits </h3>');
    fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',fullfile(outputFigures, 'Fit.html'),PosFitSpec(4)*1.5);
end

if MRSCont.flags.didQuantify
    % OspreyQuantify
    fprintf(fid,'\n<h2> Osprey Quantify</h2>');
    fprintf(fid,'\n<h3> Metabolite maps (no QC applied) </h3>');
    for ff = 1 : length(files_Quantification_Maps)
        fprintf(fid,'\n %s ', names_Quantification_Maps{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_Quantification_Maps{ff});
    end

    fprintf(fid,'\n<h3> Relative CRLB maps </h3>');
    for ff = 1 : length(files_CRLB_Maps)
        fprintf(fid,'\n %s ', names_CRLB_Maps{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_CRLB_Maps{ff});
    end

    fprintf(fid,'\n<h3> QC filter maps </h3>');
    for ff = 1 : length(files_Quantification_Maps_QC)
        fprintf(fid,'\n %s ', names_Quantification_Maps_QC{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_Quantification_Maps_QC{ff});
    end

    fprintf(fid,'\n<h3> Metabolite maps (QC filtered) </h3>');
    for ff = 1 : length(files_Quantification_Maps_QCfilt)
        fprintf(fid,'\n %s ', names_Quantification_Maps_QCfilt{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="420px" frameborder="0"></iframe>',files_Quantification_Maps_QCfilt{ff});
    end


end
if MRSCont.flags.didOverview
    % OspreyOverview
    fprintf(fid,'\n<h2> Osprey Overview</h2>');
    fprintf(fid,'\n<h3> Atlas Analysis </h3>');
    for ff = 1 : length(files_atlas)
        fprintf(fid,'\n <h4> %s </h4>', names_atlas{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',files_atlas{ff},PosGlobalConc(4)*1.5);
    end

    fprintf(fid,'\n<h3> Global Concentrations </h3>');
    for ff = 1 : length(files_GlobalConc)
        fprintf(fid,'\n <h4> %s </h4>', names_GlobalConc{ff} );
        fprintf(fid,'\n<iframe src=" %s" width="100%%" height="%ipx" frameborder="0"></iframe>',files_GlobalConc{ff},PosAtlas(4)*1.5);
    end
end
fprintf(fid,'\n</body>');
    fprintf(fid,'\n</html>');
    fclose(fid);
end
function [p] = cleanup_montages(p)
    p.layout.height = 200*2;
    p.layout.width = 840*2;
    p.layout.yaxis1.scaleanchor = 'x';
    p.layout.yaxis1.scaleratio = 1;
    p.layout.yaxis1.autorange = 'reversed';
    p.layout.xaxis1.showgrid = false;
    p.layout.yaxis1.showgrid = false;
    p.layout.xaxis1.zeroline = false;
    p.layout.yaxis1.zeroline = false;
    p.layout.xaxis1.showline = false;
    p.layout.yaxis1.showline = false;
    p.layout.xaxis1.showticklabels = false;
    p.layout.yaxis1.showticklabels = false;
    p.layout.xaxis1.domain = [0, 1];
    p.layout.yaxis1.domain = [0, 1];
    p.layout.plot_bgcolor = 'white';
    p.layout.paper_bgcolor = 'white';
    p.data{1}.colorbar.x = 1.02;
    p.data{1}.colorbar.xanchor = 'left';
    p.data{1}.colorbar.y = 0.5;
    p.data{1}.colorbar.yanchor = 'middle';
    p.data{1}.colorbar.len = 0.95;
    p.data{1}.colorbar.thickness = 20; % pixels
    p.data{1}.colorbar.thicknessmode = 'pixels';
    p.layout.margin.l = 0;
    p.layout.margin.r = 60;  % Adjust based on colorbar position
    p.layout.margin.t = 0;
    p.layout.margin.b = 0;
    p.layout.margin.pad = 0;
    plotly(p); % Update the plot
end
function [p] = cleanup_seg(p)
    p.layout.plot_bgcolor = 'black';
    p.layout.paper_bgcolor = 'black';
    p.data{1}.colorbar.tickcolor = "rgb(255,255,255)";
    p.data{1}.colorbar.tickfont.color = "rgb(255,255,255)";
    p.data{1}.colorbar.ticklabelposition = "outside";
    p.data{1}.colorbar.ticks = "outside";
    p.layout.autosize = true;
    plotly(p); % Update the plot
end
function [p] = cleanup_spectra(p)
    for dd = 1 : length(p.data)
        if strcmp(p.data{dd}.type,"scatter") && length(p.data{dd}.x) == 4
           p.data{dd}.marker.size = 5;
           p.data{dd}.marker.line.width = 0;
           p.data{dd}.mode = 'lines';
           p.data{dd}.line.color = p.data{dd}.marker.color;
           p.data{dd}.line.width = 2.5;
           p.data{dd}.line.dash  = "solid";
           tempCorner3_x = p.data{dd}.x(3);
           tempCorner3_y = p.data{dd}.y(3);
           p.data{dd}.x(3)=p.data{dd}.x(4);
           p.data{dd}.y(3)=p.data{dd}.y(4);
           p.data{dd}.x(4)=tempCorner3_x;
           p.data{dd}.y(4)=tempCorner3_y;
           p.data{dd}.x(5)=p.data{dd}.x(1);
           p.data{dd}.y(5)=p.data{dd}.y(1);
        end
    end
    plotly(p); % Update the plot
end
