function osp_onSeg( ~, ~ ,gui)
%% osp_onSeg
%   Callback function on segment button click.
%
%
%   USAGE:
%       osp_onSeg( ~, ~ ,gui);
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
%%% 1. DELETE OLD PLOTS %%%
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    set(gui.figure,'HandleVisibility','off');
    if isfield(gui.Results, 'coreg')
        if length(gui.layout.coregTab.Children) == 3
                delete( gui.layout.coregTab.Children(1) );
                delete( gui.layout.coregTab.Children(1) );
                delete( gui.layout.coregTab.Children(1) );
        else if length(gui.layout.coregTab.Children) == 2
                delete( gui.layout.coregTab.Children(1) );
                delete( gui.layout.coregTab.Children(1) );
            else
                delete( gui.layout.coregTab.Children(1) );
            end
        end
    end
    gui.layout.tabs.Selection  = 4;
    [gui,MRSCont] = osp_processingWindow(gui,MRSCont);
%%% 2. CALL OSPREYSEG %%%    
    MRSCont = OspreySeg(MRSCont);
    delete(gui.layout.dummy);
    setappdata(gui.figure,'MRSCont',MRSCont);   % Write MRSCont into hidden container in gui class 
    set(gui.figure,'HandleVisibility','on');
%%% 3. INITIALIZE OUTPUT WINDOW %%%    
    osp_iniCoregWindow(gui);       
    gui.layout.b_segm.Enable = 'off';
    if MRSCont.flags.didQuantify && ~MRSCont.flags.speedUp %Rerun quantify
        MRSCont = OspreyQuantify(MRSCont);
        gui.quant.Number.Models = length(fieldnames(MRSCont.quantify.tables));
        gui.quant.Names.Model = fieldnames(MRSCont.quantify.tables);
        gui.quant.Number.Quants = length(fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1})));
        gui.quant.Names.Quants = fieldnames(MRSCont.quantify.tables.(gui.quant.Names.Model{1}));
        gui.quant.Number.Metabs = length(MRSCont.quantify.metabs);
        gui.pop_quantOvPlotControls.String = gui.quant.Names.Quants;
        gui.pop_distrOvQuantControls.String = gui.quant.Names.Quants;
        gui.pop_corrOvQuantControls.String = gui.quant.Names.Quants;
    end
end % onSeg