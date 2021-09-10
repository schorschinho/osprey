% op_getSNR.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [SNR]=op_get_Multispectra_SNR(in,NAAppmmin,NAAppmmax,noiseppmmin,noiseppmmax);
%
% DESCRIPTION:
% Find the SNR of the NAA peak in a spectrum.
%
% INPUTS:
% in             = input data in matlab structure format
%
% OUTPUTS:
% SNR            = Estimated SNR of the input spectrum.

function [SNR]=op_get_Multispectra_SNR(in,MM);

if nargin < 2
    MM = 0;
end
% Set up noise window
noiseppmmax=0;
noiseppmmin=-2;

if ~MM
    if in.flags.isUnEdited
        SNRRange = {[1.8,2.2]};
    end
    if in.flags.isMEGA
        SNRRange = {[1.8,2.2],[2.8,3.2],[2.8,3.2],[2.8,3.2]};
    end
    if in.flags.isHERMES || in.flags.isHERCULES
        SNRRange = {[1.8,2.2],[2.8,3.2],[2.8,3.2],[1.8,2.2],[2.8,3.2],[2.8,3.2],[2.8,3.2]};
    end
else
    if in.flags.isUnEdited
        SNRRange = {[0.7,1.1]};
    end
    if in.flags.isMEGA
        SNRRange = {[0.7,1.1],[0.7,1.1],[0.7,1.1],[0.7,1.1]};
    end
    if in.flags.isHERMES || in.flags.isHERCULES
        SNRRange = {[0.7,1.1],[0.7,1.1],[0.7,1.1],[0.7,1.1],[0.7,1.1],[0.7,1.1],[0.7,1.1]};
    end   
end

SNR = cell(in.subspecs,1);
for ss = 1 : in.subspecs
    if in.subspecs == 1
        temp{ss} = in;
    else
        temp{ss} = op_takesubspec(in,ss);
    end
    NAAppmmin = SNRRange{ss}(1);
    NAAppmmax = SNRRange{ss}(2);
    
    %FIRST FIND THE NAA SIGNAL INTENSITY.  USE THE MAX PEAK HEIGHT OF THE
    %MAGNITUDE SPECTRUM INSIDE THE DESIRED SPECTRAL RANGE:
    NAAwindow=temp{ss}.specs(temp{ss}.ppm>NAAppmmin & temp{ss}.ppm<NAAppmmax);
    ppmwindow=temp{ss}.ppm(temp{ss}.ppm>NAAppmmin & temp{ss}.ppm<NAAppmmax);

    maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
    maxNAA=abs(NAAwindow(maxNAA_index));

    %figure;
    %plot(ppmwindow,abs(real(NAAwindow)));

    %figure
    %plot(in.ppm,real(in.specs));
    %noiseppmmin=input('input lower ppm limit for noise: ');
    %noiseppmmax=input('input upper ppm limit for noise: ');

    %NOW FIND THE STANDARD DEVIATION OF THE NOISE:
    noisewindow=temp{ss}.specs(temp{ss}.ppm>noiseppmmin & temp{ss}.ppm<noiseppmmax);
    ppmwindow2=temp{ss}.ppm(temp{ss}.ppm>noiseppmmin & temp{ss}.ppm<noiseppmmax)';

    P=polyfit(ppmwindow2,noisewindow,2);
    noise=noisewindow-polyval(P,ppmwindow2);
    %figure
    %plot(ppmwindow2,real(noisewindow),...
    %    ppmwindow2,real(polyval(P,ppmwindow2)),...
    %    ppmwindow2,real(noise));

    signal=(maxNAA-mean(real(noisewindow))); %Removes DC offset

    noisesd=std(real(noise));

    %SNR=maxNAA/noisesd
    SNR{ss}=signal/noisesd;

    %For out of volume MRSI
    if isempty(SNR{ss}) || isnan(sum(SNR{ss})) || length(SNR{ss}) > 1
        SNR{ss} = 0;
    end
end
