function MRSCont = osp_interactive_mrsi_mask(MRSCont,kk,yLim)
if nargin < 2
    kk =1;
end
if nargin < 3
    yLim = [-0.0001 0.001]; % Set y limits
end
    Norm = MRSCont.plot.processed.w.ContMax;

    if MRSCont.processed.w{kk}.nZvoxels > 1
        mask = zeros(MRSCont.processed.w{kk}.nXvoxels,MRSCont.processed.w{kk}.nYvoxels,MRSCont.processed.w{kk}.nZvoxels);
    else
        mask = zeros(MRSCont.processed.w{kk}.nXvoxels,MRSCont.processed.w{kk}.nYvoxels);
    end

    for z = 1 : MRSCont.processed.w{kk}.nZvoxels
        scales = squeeze(max(abs(real(MRSCont.processed.w{kk}.specs(:,:,:,z)))/Norm));
        
        scales(scales > 0.07) = 1;
        scales(scales < 0.07) = 0;
        mask(:,:,z) = scales;
        ppmmin = 0;
        ppmmax = 4;
%         ppmmin = 4.65-1;
%         ppmmax = 4.65+1;
        add_No_MoCo = 0;% Include spectra without MoCo in the plot?
        coord = [1 MRSCont.processed.w{kk}.nXvoxels; 1 MRSCont.processed.w{kk}.nYvoxels]; % Which coordinates
        
        which_spec = {'A'};
%         out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, 0,0,z,yLim,ones(coord(2),coord(1)),1);
        out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, 0,0,z,yLim,mask(:,:,z),0);
        size = get(gcf, 'Position');
        set(gcf,'Position',[size(1) size(2) 1.2*size(3) 1.2*size(3)*1.2])
        set(gcf,'Renderer','painters','Menu','none','Toolbar','none');
        mkdir(fullfile(MRSCont.outputFolder,'Figures'));
        saveas(out,fullfile(MRSCont.outputFolder,'Figures',['slice_' num2str(z) '_water']),'pdf');
        close(gcf);
        open(fullfile(MRSCont.outputFolder,'Figures',['slice_' num2str(z) '_water.pdf']));
%         out3 = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, 0,0,z,yLim,mask(:,:,z),0);
%         size = get(gcf, 'Position');
%         set(gcf,'Position',[size(1) size(2) 1.2*size(3) 1.2*size(3)*1.2])
%         set(gcf,'Renderer','painters','Menu','none','Toolbar','none');
        
        out2 = figure;
       imagesc(mask(:,:,z)');
       title('Move ellipse to cover voxels of interst. Double click when finished');
       size = get(gcf, 'Position');
       set(gca, 'XTick',1:MRSCont.processed.w{kk}.nXvoxels,'YTick',1:MRSCont.processed.w{kk}.nYvoxels,'TickDir','out','XGrid','on','YGrid','on','GridAlpha',1)
        set(gcf,'Position',[size(1) size(2) 1.2*size(3) 1.2*size(3)*1.2])
        set(gcf,'Renderer','painters','Menu','none','Toolbar','none');
       h = drawellipse('Center',[8 9],'SemiAxes',[6 8], ...
    'StripeColor','y');
       pos = customWait(h);
       mask(:,:,z) = createMask(h)';
       close(gcf)
        out2 = figure;
       imagesc(mask(:,:,z)');
       title('Add/remove voxels by clicking. Press return when finished.');
       size = get(gcf, 'Position');
       set(gca, 'XTick',1:MRSCont.processed.w{kk}.nXvoxels,'YTick',1:MRSCont.processed.w{kk}.nYvoxels,'TickDir','out','XGrid','on','YGrid','on','GridAlpha',1)
        set(gcf,'Position',[size(1) size(2) 1.2*size(3) 1.2*size(3)*1.2])
        set(gcf,'Renderer','painters','Menu','none','Toolbar','none');
        while true
             x= [];
            y =[];
           [x, y, button] = ginput(1);
           if isempty(x); break; end
            x= round(x);
           y = round(y);
           for i = 1 : length(x)
                if mask(x(i),y(i),z) == 1
                   mask(x(i),y(i),z) = 0;
                else
                   mask(x(i),y(i),z) = 1;
                end
           end
           imagesc(mask(:,:,z)');
           drawnow
           set(gca, 'XTick',1:MRSCont.processed.w{kk}.nXvoxels,'YTick',1:MRSCont.processed.w{kk}.nYvoxels,'TickDir','out','XGrid','on','YGrid','on','GridAlpha',1)
        end
        close(gcf);

       if z == 2
           out2 = figure;
           imagesc(mask(:,:,z)');
           title('Add center voxel for spiral by clicking. Press return when finished.');
           size = get(gcf, 'Position');
       set(gca, 'XTick',1:MRSCont.processed.w{kk}.nXvoxels,'YTick',1:MRSCont.processed.w{kk}.nYvoxels,'TickDir','out','XGrid','on','YGrid','on','GridAlpha',1)
        set(gcf,'Position',[size(1) size(2) 1.2*size(3) 1.2*size(3)*1.2])
        set(gcf,'Renderer','painters','Menu','none','Toolbar','none');
        while true
           x= [];
           y =[];
           [x,y] = ginput(1);
           if isempty(x); break; end
           x= round(x);
           y = round(y);
           for i = 1 : length(x)
                   mask(x(i),y(i),z) = 2;
           end
           imagesc(mask(:,:,z)');
           drawnow
           set(gca, 'XTick',1:MRSCont.processed.w{kk}.nXvoxels,'YTick',1:MRSCont.processed.w{kk}.nYvoxels,'TickDir','out','XGrid','on','YGrid','on','GridAlpha',1)
        end
           close(gcf);
       end
    end
   
    MRSCont.mask{kk} = mask;

    outputFolder = MRSCont.outputFolder;
    outputFile      = MRSCont.outputFile;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end