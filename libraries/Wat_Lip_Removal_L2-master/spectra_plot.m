function [] = spectra_plot(x_start, y_start, x_size, y_size, p_first, p_end, line, csi1, csi2, csi3, csi4)
% spectra plot for csi data
% csi : input spectra
% x_start, y_start : upper left corner of ROI
% x_size, y_size : size of ROI to show
% p_first, p_end : indices of spectral region by ppm, optional,
% line: show peak lines, Y(1)/N(0)

N = size(csi1);

% transfer from ppm to data points
f_first = round((12.5-p_end)/15.6*N(3));
f_end = round((12.5-p_first)/15.6*N(3));

% line position, optional
s1=round(0.5903*N(3));
e1=round(0.6028*N(3));
s2=round(0.6028*N(3));
e2=round(0.6153*N(3));
s3=round(0.6670*N(3));
e3=round(0.6795*N(3));
s4=round(0.5453*N(3));
e4=round(0.5578*N(3));
s5=round(0.6460*N(3));
e5=round(0.6585*N(3));
startp=[s1 s2 s3 s4 s5];
endp=[e1 e2 e3 e4 e5];

% set the maximum in y axis
csi_ = abs(csi1(x_start:x_start+x_size-1, y_start:y_start+y_size-1, f_first:f_end));
y_max = 1.5*max(csi_(:));
scale = 0.25*y_max;

% figure layout
nnx = x_size;
nny = y_size;
widthx = 1/nnx*.9;
widthy = 1/nny*.9;
marginx = 1/nnx*0.05;
marginy = 1/nny*0.05;
originx = marginx:1/nnx:1;
originy = 1-1/nny+marginy:-1/nny:0;

xorig = x_start;
yorig = y_start;
x = xorig:xorig+nnx-1;
y = yorig:yorig+nny-1;

figure();
for j=1:nny
    for i=1:nnx
        subplot('Position', [originy(nny-j+1), originx(nnx-i+1), widthy, widthx])
        hold on, plot(f_first:f_end, squeeze(csi1(x(i),y(j),f_first:f_end)),'b', 'LineWidth', 1);
        if nargin == 9
            hold on, plot(f_first:f_end, squeeze(csi2(x(i),y(j),f_first:f_end))+scale,'r', 'LineWidth', 1);
        end
        if nargin == 10
            hold on, plot(f_first:f_end, squeeze(csi2(x(i),y(j),f_first:f_end))+scale,'r', 'LineWidth', 1);
            hold on, plot(f_first:f_end, squeeze(csi3(x(i),y(j),f_first:f_end))+2*scale,'g', 'LineWidth', 1);
        end
        if nargin == 11
            hold on, plot(f_first:f_end, squeeze(csi2(x(i),y(j),f_first:f_end))+scale,'r', 'LineWidth', 1);
            hold on, plot(f_first:f_end, squeeze(csi3(x(i),y(j),f_first:f_end))+2*scale,'g', 'LineWidth', 1);
            hold on, plot(f_first:f_end, squeeze(csi4(x(i),y(j),f_first:f_end))+3*scale,'k', 'LineWidth', 1);
        end
        
        % draw reference lines
        if line == 1
%             maxval=max(abs(squeeze(csi1(x(i),y(j), f_first:f_end))));
            for kk=1:4
                hold on, plot(ones(1,100)*startp(kk),linspace(0,y_max/1,100),'r','LineWidth', 0.5);
                hold on, plot(ones(1,100)*endp(kk),linspace(0,y_max/1,100),'r','LineWidth', 0.5);
            end
        end
%         axis([f_first,f_end, -scale/5, y_max/1]);
        axis([1000, 1500, -60 300])
    end
end
hold off;
end


