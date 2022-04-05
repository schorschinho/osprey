function MRSCont = osp_interactive_mrsi_mask(MRSCont,kk)
if nargin < 2
    kk =1;
end

    spec = op_freqrange(MRSCont.raw{1},4.0,5.5);
    mask = squeeze(sum(squeeze(abs(real(spec.specs))),1));
    mask = mask/(max(max(max(mask))));
    mask(mask > 2 * mask(1,1)) = 1;
    mask(mask < 1) = 0;  
    
   out = figure('WindowButtonDownFcn',@callBack);
   heatmap(mask');
   
%    h = drawellipse;
%    roi = images.roi.Freehand(gca,'Position',[1 1;1 4;4 1;4 4]);

end
 function callBack(hObject,~)
    mousePos=get(hObject,'CurrentPoint');
    disp(['You clicked X:',num2str(mousePos(1)),',  Y:',num2str(mousePos(2))]);
 end