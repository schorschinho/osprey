function osp_updatemeanOvWindow(gui)
%% osp_updatemeanOvWindow
%   This function updates the mean spectra overview tab.
%
%
%   USAGE:
%       osp_updatemeanOvWindow(gui);
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
        delete(gui.Plot.meanOv.Children(2).Children)
%%% 2. VISUALIZATION PART OF THIS TAB %%%
        if (MRSCont.flags.isMEGA) %Is Edited? Pick the right tab
            if (gui.process.Selected == 1 || gui.process.Selected == 2 || gui.process.Selected == 3)
                t = gui.process.Selected;
                if (strcmp(gui.quant.Names.Model,'conc'))
                    fit = 1;
                end                    
            else if gui.process.Selected == 4
                    t = 5;
                    if (strcmp(gui.quant.Names.Model,'conc'))
                        fit = 2;
                    end                     
                else if MRSCont.flags.hasWater
                    t = 6;
                    if (strcmp(gui.quant.Names.Model,'conc'))
                        fit = 3;
                    end 
                    end
                end
            end
        else
            t = gui.process.Selected;
        end
        for g = 1 :  gui.overview.Number.Groups
            temp = figure;
            hold on
            if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
                %Metabolite data
                if (strcmp(gui.quant.Names.Model,'conc') && ~strcmp(gui.process.Names{t},'A'))%Is Concatenated 
                    if gui.overview.NAAnormed ==1 %Is NAA normalized
                        shift = gui.layout.shift * (g-1);
                        yu = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) + ...
                            MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                        yl = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) - ...
                            MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                        temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_conc_fit' num2str(fit)])+shift ,'color', MRSCont.colormap.Accent, 'LineWidth', 1);
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_conc_baseline_' num2str(fit)])+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1);
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_conc_res_' num2str(fit)])+shift ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);
                        plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 2);                        
                    else %Not NAA normalized ... Normalize to max
                        ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.process.Names{2}]));
                        shift = ylimmax * gui.layout.shift * (g-1);
                        yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) + ...
                            MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                        yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) - ...
                            MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                        temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.Accent, 'LineWidth', 1);
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_conc_baseline_' num2str(fit)])+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1);
                        plot(MRSCont.overview.ppm_fit_conc,MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_conc_res_' num2str(fit)])+shift ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);
                        plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_conc_fit' num2str(fit)])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 1);
                    end                
                else
                    if gui.overview.NAAnormed ==1 %Is NAA normalized
                        shift = gui.layout.shift * (g-1);
                        yu = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) + ...
                            MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                        yl = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) - ...
                            MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                        temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                        if (strcmp(gui.quant.Names.Model,'conc') && strcmp(gui.process.Names{t},'A'))%Is Concatenated 
                            plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 2);
                        else
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.Accent, 'LineWidth', 1);
                            plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 2);
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_baseline_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1);
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g) '_NAAnormalized']).(['mean_res_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);
                        end
                    else %Not NAA normalized ... Normalize to max
                        ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.process.Names{2}]));
                        shift = ylimmax * gui.layout.shift * (g-1);
                        yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) + ...
                            MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                        yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) - ...
                            MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                        temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                        if (strcmp(gui.quant.Names.Model,'conc') && strcmp(gui.process.Names{t},'A'))%Is Concatenated 
                            plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 1);
                        else
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.Accent, 'LineWidth', 1);
                            plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 1);
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_baseline_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1);
                            plot(MRSCont.overview.(['ppm_fit_' (gui.fit.Names{t})]),MRSCont.overview.(['sort_fit_g' num2str(g)]).(['mean_res_' gui.fit.Names{t}])+shift ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);
                        end
                    end
                end 
            else if (strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D'))
                if gui.overview.NAAnormed ==1 %Is NAA normalized
                    shift = gui.layout.shift * (g-1);
                    yu = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) + ...
                        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                    yl = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}]) - ...
                        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.process.Names{t}]);
                    temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                    plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 2);
                else %Not NAA normalized ... Normalize to max
                    ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.process.Names{2}]));
                    shift = ylimmax * gui.layout.shift * (g-1);
                    yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) + ...
                        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                    yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) - ...
                        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                    temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                    plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 2);
                end
            else %Water data
                ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.process.Names{1}]));
                shift = ylimmax * gui.layout.shift * (g-1);
                yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) + ...
                    MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}]) - ...
                    MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.process.Names{t}]);
                temp = fill([MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]) fliplr(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                plot(MRSCont.overview.(['ppm_data_' (gui.process.Names{t})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.process.Names{t}])+shift ,'color',gui.colormap.cb(g,:), 'LineWidth', 1);
                end
            end
                figpl = get(temp.Parent,'Parent');
                set(temp.Parent.Children,'Parent',gui.Plot.meanOv.Children(2))
                close(figpl);
        end
        if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
            set(gui.Plot.meanOv.Children(2), 'XLim', [0.2 4.5])
            set(gui.Plot.meanOv.Children(2).Title, 'String', ['Mean \pm SD ' gui.layout.proTab.TabTitles{gui.process.Selected}])
        else
            set(gui.Plot.meanOv.Children(2), 'XLim', [0 2*4.68])
            set(gui.Plot.meanOv.Children(2).Title, 'String', ['Mean \pm SD ' gui.layout.proTab.TabTitles{gui.process.Selected}])
        end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class        
end