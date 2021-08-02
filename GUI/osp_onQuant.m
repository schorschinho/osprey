function osp_onQuant( ~, ~ ,gui)
%% osp_onQuant
%   Callback function on quantify button click.
%
%
%   USAGE:
%       osp_onQuant( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-16: First version of the code.
%%% 1. INITIALIZE %%%
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    set(gui.figure,'HandleVisibility','off');
    gui.layout.tabs.Selection  = 5;
    [gui,MRSCont] = osp_processingWindow(gui,MRSCont);
%%% 2. CALL OSPREYQUANTIFY %%%    
    MRSCont = OspreyQuantify(MRSCont); 
    gui.quant.Number.Model = length(fieldnames(MRSCont.quantify.tables));
    gui.quant.Names.Model = fieldnames(MRSCont.quantify.tables);
    gui.quant.Number.Quants = length(fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1})));
    gui.quant.Names.Quants = fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1}));
    gui.quant.Number.Metabs = length(MRSCont.quantify.metabs);
    gui.overview.Selected.Metab = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{1}), 'tNAA'));
    if isempty(gui.overview.Selected.Metab)
        gui.overview.Selected.Metab = find(strcmp(MRSCont.quantify.metabs.(gui.quant.Names.Model{1}), 'NAA'));
    end
    gui.quant.idx.GABA = find(strcmp(MRSCont.quantify.metabs, 'GABA'));
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
%%% 3. INITIALIZE QUANTIFY OUTPUT WINDOW %%% 
    MRSCont = OspreyOverview(MRSCont);
    delete(gui.layout.dummy); 
    osp_iniQuantifyWindow(gui);
    gui.overview.NAAnormed = 1;
    gui.overview.Number.Groups = MRSCont.overview.NoGroups;
    [gui.colormap.cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
    temp = gui.colormap.cb(3,:);
    gui.colormap.cb(3,:) = gui.colormap.cb(4,:);
    gui.colormap.cb(4,:) = temp;
    if isfield(MRSCont.overview, 'corr')
        gui.overview.Names.Corr = MRSCont.overview.corr.Names{1};
        gui.overview.CorrMeas = MRSCont.overview.corr.Meas;
    end
    gui.overview.Selected.Corr = 1;
    gui.overview.Selected.CorrChoice = 1;
    gui.overview.Names.QM = {'SNR','FWHM (ppm)'};
    set(gui.layout.quantifyTab, 'SelectionChangedFcn',{@osp_QuantTabChangeFcn,gui});      
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
%%% 4. CALL OSPREYOVERVIEW %%%    
    osp_iniOverviewWindow(gui);
    set(gui.layout.overviewTab, 'SelectionChangedFcn',{@osp_OverviewTabChangedFcn,gui});
    set(gui.controls.pop_specsOvPlot,'callback',{@osp_pop_specsOvPlot_Call,gui});
    set(gui.controls.pop_meanOvPlot,'callback',{@osp_pop_meanOvPlot_Call,gui});
    set(gui.controls.pop_quantOvPlot,'callback',{@osp_pop_quantOvPlot_Call,gui});
    set(gui.controls.pop_distrOvQuant,'callback',{@osp_pop_distrOvQuant_Call,gui});
    set(gui.controls.pop_distrOvMetab,'callback',{@osp_pop_distrOvMetab_Call,gui});
    set(gui.controls.pop_corrOvQuant,'callback',{@osp_pop_corrOvQuant_Call,gui});
    set(gui.controls.pop_corrOvMetab,'callback',{@osp_pop_corrOvMetab_Call,gui});
    set(gui.controls.pop_corrOvCorr,'callback',{@osp_pop_corrOvCorr_Call,gui});
    set(gui.controls.pop_whichcorrOvCorr,'callback',{@osp_pop_whichcorrOvCorr_Call,gui});
    gui.layout.b_quant.Enable = 'off';
    set(gui.figure,'HandleVisibility','on');
end % onQuant