function wat = generator_water(SW, Datapoints, Zerofill, wat_frq_range, damp, damp_range, Watnum, echoposition)
%% Parameters
%SW = SW;  % spectral width
NF2 = Datapoints; % data points
NNW = Watnum; % number of water components
WR = wat_frq_range;   % water frequency dispersion range
t = 1/SW:1/SW:1/SW*NF2; % time points

%% generate fid
fid = zeros(NF2,NNW);
for ii = 1:NNW
    %% frequence, damping facter, phase, and amplitude
    feq_wat = 2*WR*(rand(1)-0.5);
    dam_wat = (rand(1))*(damp_range-damp)+damp;
    pha_wat = (rand(1)-0.5)*2*pi;
    amp_wat = 0.5; %%*exp(-feq_wat/WR);
    %% signal model
        fid(:,ii) = amp_wat*exp(-dam_wat*abs(t-echoposition*NF2*1/SW)).*exp(-i*((feq_wat)*2*pi*t-pha_wat));
end
wat = flipud(fftshift((fft(fid,Zerofill))));
%% plot figures
% figure();
% subplot(2,1,1); plot(real(fid),'linewidth',1)
% subplot(2,1,2); plot(real(wat),'linewidth',1)
end