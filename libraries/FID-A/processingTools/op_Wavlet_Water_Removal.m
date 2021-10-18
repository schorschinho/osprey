function out = op_Wavlet_Water_Removal(in)
% Based on: Cobas, J.C., Bernstein, M.A., Martín-Pastor, M., Tahoces, P.G.,
% 2006. A new general-purpose fully automatic baseline-correction procedure
% for 1D and 2D NMR data. J. Magn. Reson. 183, 145-51.

out = in;

showPlots =0;

[z.real, z.imag] = BaselineModeling(in);

if showPlots == 1
    
    freq = in.ppm;
    
    H = figure;
    d.w = 0.5;
    d.h = 0.65;
    d.l = (1-d.w)/2;
    d.b = (1-d.h)/2;
    set(H,'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [d.l d.b d.w d.h]);
    
    subplot(2,2,1); cla;
    plot(freq, real(in.specs), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    hold on;
    plot(freq, z.real, 'r');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,2); cla;
    plot(freq, imag(in.specs), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    hold on;
    plot(freq, z.imag, 'r');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,3); cla;
    plot(freq, real(in.specs) - z.real, 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,4); cla;
    plot(freq, imag(in.specs) - z.imag, 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    drawnow;
    %pause(0.5);
    
end

out.specs = complex(real(in.specs) - z.real, imag(in.specs) - z.imag);
out.fids = ifft(fftshift(complex(real(in.specs) - z.real, imag(in.specs) - z.imag),in.dims.t),[],in.dims.t);

out.watersupp.fids_water = complex(z.real, z.imag);
out.watersupp.specs_water = ifft(fftshift(complex(z.real,z.imag),in.dims.t),[],in.dims.t);
out.watersupp.amp = max(real(out.watersupp.specs_water));


end


function [z_real, z_imag] = BaselineModeling(in)
% Based on:
% Golotvin & Williams, 2000. Improved baseline recognition
%   and modeling of FT NMR spectra. J. Magn. Reson. 146, 122-125
% Cobas et al., 2006. A new general-purpose fully automatic
%   baseline-correction procedure for 1D and 2D NMR data. J. Magn. Reson.
%   183, 145-151

% Power spectrum of first-derivative of signal calculated by CWT

spec = in.specs;

Wy = abs(cwt2(real(spec), 10)).^2;

freq = in.ppm;
noiseLim = freq <= 0 & freq >= -2;

sigma = std(Wy(noiseLim));

w = 1:5;
k = 3;
baseline = zeros(length(Wy),1);
signal = zeros(length(Wy),1);

while 1
    if w(end) > length(Wy)
        break
    end
    h = max(Wy(w)) - min(Wy(w));
    if h < k*sigma
        baseline(w) = Wy(w);
    else
        signal(w) = Wy(w);
    end
    w = w + 1;
end

% Include lipids in baseline estimate and water, as appropriate

waterLim = freq <= 8.5 & freq >= 4.25;
baseline(waterLim) = Wy(waterLim);

z_real = real(spec);
z_real(baseline == 0) = 0;

water = whittaker(z_real(waterLim), 2, 0.2);
z_real = whittaker(z_real, 2, 1e3);

z_real(waterLim) = water;

z_imag = -imag(hilbert(z_real));

end


function Wy = cwt2(y, a)

precis = 10;
coef = sqrt(2)^precis;
pas  = 1/2^precis;

lo    = [sqrt(2)*0.5 sqrt(2)*0.5];
hi    = lo .* [1 -1];
long  = length(lo);
nbpts = (long-1)/pas+2;

psi = coef*upcoef2(lo,hi,precis);
psi = [0 psi 0];
xVal = linspace(0,(nbpts-1)*pas,nbpts);

step = xVal(2) - xVal(1);
valWAV = cumsum(psi)*step;
xVal = xVal - xVal(1);
xMaxVal = xVal(end);

y = y(:)';
j = 1+floor((0:a*xMaxVal)/(a*step));
if length(j) == 1
    j = [1 1];
end
f = fliplr(valWAV(j));

Wy = -sqrt(a) * keepVec(diff(conv2(y,f,'full')),length(y));

    function out = upcoef2(lo,hi,a)
        out = hi;
        for k = 2:a
            out = conv2(dyadup2(out),lo,'full');
        end
        function out = dyadup2(in)
            z = zeros(1,length(in));
            out = [in; z];
            out(end) = [];
            out = out(:)';
        end
    end

    function out = keepVec(in,len)
        out = in;
        ok = len >= 0 && len < length(in);
        if ~ok
            return
        end
        d     = (length(in)-len)/2;
        first = 1+floor(d);
        last  = length(in)-ceil(d);
        out   = in(first:last);
    end

end

function z = whittaker(y, d, lambda)
% Code taken from: Eilers, 2003. A perfect smoother. Anal. Chem. 75,
% 3631-3636

y = y(:);

m = length(y);
E = speye(m);
D = diff(E,d);

w = double(y ~= 0);
W = spdiags(w,0,m,m);
C = chol(W + lambda*(D'*D));
z = C\(C'\(w.*y));

end