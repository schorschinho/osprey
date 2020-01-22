function osp_updateSpecsOvWindow(gui)
%% osp_updateSpecsOvWindow
%   This function updates the specs overview tab.
%
%
%   USAGE:
%       osp_updateSpecsOvWindow(gui);
%
%   INPUT:  
%           gui      = gui class containing all handles and the MRSCont             
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
        MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class
        delete(gui.Plot.specsOv.Children(2).Children)
%%% 2. VISUALIZATION PART OF THIS TAB %%%
        if (MRSCont.flags.isMEGA) %Is Edited? Pick the right tab
            if (gui.process.Selected == 1 || gui.process.Selected == 2 || gui.process.Selected == 3)
                t = gui.process.Selected;
            else if gui.process.Selected == 4
                    t = 5; 
                    else if gui.process.Selected == 5
                    t = 4;
                    else if MRSCont.flags.hasWater
                        t = 6;
                        end
                    end
                end
            end
        else
            t = gui.process.Selected;
        end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is Edited? Pick the right tab
            if (gui.process.Selected == 1 || gui.process.Selected == 2 || gui.process.Selected == 3 || gui.process.Selected == 4 || gui.process.Selected == 5 || gui.process.Selected == 6)
                t = gui.process.Selected;
            else if gui.process.Selected == 7
                    t = 8; 
                    else if gui.process.Selected == 8
                    t = 7;
                    else if MRSCont.flags.hasWater
                        t = 9;
                        end
                    end
                end
            end
        end
        for g = 1 :  gui.overview.Number.Groups %Loop over groups
            temp = figure( 'Visible', 'off' );
            if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
                ylimmax = max(real(MRSCont.overview.all_data.(gui.process.Names{t}){1,1}.specs));
                shift = ylimmax * gui.layout.shiftind * (g-1);
                temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.process.Names{t}),2,1,gui.colormap.cb(g,:),shift,['Overview ' gui.layout.proTab.TabTitles{t}]);
            else %Water data?
                ylimmax = max(real(MRSCont.overview.all_data.(gui.process.Names{1}){1,1}.specs));
                shift = ylimmax * gui.layout.shiftind * (g-1);
                temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.process.Names{t}),2,1,gui.colormap.cb(g,:),shift,['Overview ' gui.layout.proTab.TabTitles{t}]);
            end
            set(gca, 'YColor', MRSCont.colormap.Background);
            set(gca,'YTickLabel',{})
            set(gca,'YTick',{});
            set(gca,'XColor',MRSCont.colormap.Foreground);
            set(gca,'Color','w');
            set(gcf,'Color','w');
            title(['Overview ' gui.layout.proTab.TabTitles{t}],'Color', MRSCont.colormap.Foreground);
                ax=get(temp,'Parent');
                figpl = get(ax,'Parent');
                copyobj(ax.Children, gui.Plot.specsOv.Children(2));
                close( figpl );
        end
        if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
            set(gui.Plot.specsOv.Children(2), 'XLim', [0.2 4.5])
            set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' gui.layout.proTab.TabTitles{gui.process.Selected}])
        else
            set(gui.Plot.specsOv.Children(2), 'XLim', [0 2*4.68])
            set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' gui.layout.proTab.TabTitles{gui.process.Selected}])
        end
end
