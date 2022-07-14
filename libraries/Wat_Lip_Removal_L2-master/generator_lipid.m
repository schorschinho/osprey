function lip = generator_lipid(SW, Datapoints, Zerofill, Lipid_Left, Lipid_Right, damp, damp_range, spe_number, echopos)
%% Parameters
NF2 = Datapoints; % data points
NNL = spe_number; % lipid components
LR = Lipid_Right;   % lipid frequencey dispersion
LF = Lipid_Left;   %lipid frequency
t = 1/SW:1/SW:NF2/SW; % time points

%% generate fid
fid = zeros(NF2,NNL);
for ii = 1:NNL
    %% frequence, damping facter, phase, and amplitude
    feq_lip = LF+(LR-LF)*rand(1);
    dam_lip = (rand(1))*(damp_range-damp)+damp;
    pha_lip = (rand(1)-0.5)*2*pi;
    amp_lip = exp(-((abs(feq_lip-(LR+LF)/2)).^2)/80^2);
%     amp_lip = 1;
%     amp_lip = 100;
    %% signal model
    fid(:,ii) = amp_lip*exp(-dam_lip*abs(t-echopos*NF2*1/SW)).*exp(-i*((feq_lip)*2*pi*t-pha_lip));
end
lip = flipud(fftshift((fft(fid,Zerofill))));
end