function MM_spectrum = MM_spectrum_model(input,ppm)
%ppm is the axis that defines the ppm sampling of MM_spectrum
%Spectrum is based upon shifts in Giapitzakis et al. 

%input(1) is an amplitude and input(2) a width, etc.
%input29 is a slope and input30 an offset
%last input parameter input(end) is a ppm referencing shift 



ppm_09=0.94+input(end);
ppm_12=1.22+input(end);
ppm_14=1.43+input(end);
ppm_17=1.69+input(end);
ppm_20=2.04+input(end);
ppm_22=2.27+input(end);
ppm_27=2.66+input(end);%change from Giap due to stronger coupling at 3T
ppm_30=3.01+input(end);
ppm_32=3.21+input(end);
ppm_37=3.71+input(end);
ppm_38=3.79+input(end);
ppm_39=3.87+input(end);
ppm_40=3.97+input(end);
ppm_42=4.2+input(end);

MM_spectrum = GaussModel([input(1) input(2) ppm_09 0 0],ppm)+GaussModel([input(3) input(4) ppm_12 0 0],ppm)+GaussModel([input(5) input(6) ppm_14 0 0],ppm)+GaussModel([input(7) input(8) ppm_17 0 0],ppm)+GaussModel([input(9) input(10) ppm_20 0 0],ppm)+GaussModel([input(11) input(12) ppm_22 0 0],ppm)+GaussModel([input(13) input(14) ppm_27 0 0],ppm)+GaussModel([input(15) input(16) ppm_30 0 0],ppm)+GaussModel([input(17) input(18) ppm_32 0 0],ppm)+GaussModel([input(19) input(20) ppm_37 0 0],ppm)+GaussModel([input(21) input(22) ppm_38 0 0],ppm)+GaussModel([input(23) input(24) ppm_39 0 0],ppm)+GaussModel([input(25) input(26) ppm_40 0 0],ppm)+GaussModel([input(27) input(28) ppm_42 0 0],ppm)+input(29)*(ppm-2.5)+input(30);

% MM_spectrum = MM_spectrum.';



end