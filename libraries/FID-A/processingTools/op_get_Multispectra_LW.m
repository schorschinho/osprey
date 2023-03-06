% op_get_Multispectra_LW.m
% Helge Zoellner, Johns Hopkins University 2014.
%
% USAGE:
% [FWHM]=op_get_Multispectra_LW(in);
%
% DESCRIPTION:
% Estimates the linewidth of a reference peak in the spectrum.  By default,
% the reference peak is water, between 4.4 and 5.0 ppm.  Two methods are
% used to estimate the linewidth:  1.  FWHM is measured by simply taking the
% full width at half max of the reference peak.  2.  The FWHM is measured by
% fitting the reference peak to a lorentzian lineshape and determine the FWHM of the
% best fit.  The output FWHM is given by the average of these two measures.
%
% INPUTS:
% in         = input spectrum in structure format.
%
% OUTPUTS:
% FWHM       = Estimated linewidth of the input spectrum (in Hz).


function [FWHM]=op_get_Multispectra_LW(in);
if nargin < 2
    MM = 0;
end

zpfactor=8;
if ~MM
    if in.flags.isUnEdited
        SNRRange = {[2.9,3.1]};
    end
    if in.flags.isMEGA
        SNRRange = {[2.9,3.1],[2.9,3.1],[2.8,3.2],[2.9,3.1]};
    end
    if in.flags.isHERMES || in.flags.isHERCULES
        if size(in.target,1) == 2
            SNRRange = {[2.9,3.1],[2.9,3.1],[2.9,3.1],[2.9,3.1],[2.8,3.2],[2.8,3.2],[2.9,3.1]};
        else
            SNRRange = {[2.9,3.1],[2.9,3.1],[2.9,3.1],[2.9,3.1],[2.8,3.2],[2.8,3.2],[0.9,1.5],[2.9,3.1]};
        end
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


FWHM = zeros(in.subspecs,1);
for ss = 1 : in.subspecs
        if in.subspecs == 1
            temp{ss} = in;
        else
            temp{ss} = op_takesubspec(in,ss);
        end
    Refppmmin = SNRRange{ss}(1);
    Refppmmax = SNRRange{ss}(2);
    
    temp{ss}=op_zeropad(temp{ss},zpfactor);

    %FIRST FIND FWHM USING TWO METHODS:

    %METHOD 1:  ACTUALLY MEAUSURE FWHM OF WATER PEAK
    Refwindow=temp{ss}.specs(temp{ss}.ppm>Refppmmin & temp{ss}.ppm<Refppmmax);
    ppmwindow=temp{ss}.ppm(temp{ss}.ppm>Refppmmin & temp{ss}.ppm<Refppmmax);

    maxRef_index=find(abs(real(Refwindow))==max(abs(real((Refwindow)))));
    maxRef=real(Refwindow(maxRef_index));

    % plot(ppmwindow,abs(real(Refwindow)),'.');
    if ~isempty(maxRef) && ~(sum(maxRef) == 0) && ~(length(maxRef) > 1)
        gtHalfMax=find(abs(real(Refwindow)) >= 0.5*abs(maxRef));


        FWHM1=abs(ppmwindow(gtHalfMax(1)) - ppmwindow(gtHalfMax(end)));
        FWHM1=FWHM1*(42.577*temp{ss}.Bo);  %Assumes proton.


        %METHOD 2:  FIT WATER PEAK TO DETERMINE FWHM PARAM
        sat='n';
        waterFreq=ppmwindow(maxRef_index);
        while sat=='n'
            parsGuess=zeros(1,5);
            parsGuess(1)=maxRef; %AMPLITUDE
            parsGuess(2)=(5*temp{ss}.Bo/3)/(42.577*temp{ss}.Bo); %FWHM.  Assumes Proton.  LW = 5/3 Hz/T.
            parsGuess(3)=waterFreq; %FREQUENCY
            parsGuess(4)=0; %Baseline Offset
            parsGuess(5)=0; %Phase

            yGuess=op_lorentz(parsGuess,ppmwindow);
            parsFit=nlinfit(ppmwindow,real(Refwindow'),@op_lorentz,parsGuess);
            yFit=op_lorentz(parsFit,ppmwindow);

        %     figure;
        %     plot(ppmwindow,Refwindow,'.',ppmwindow,yGuess,':',ppmwindow,yFit);
        %     legend('data','guess','fit');
        %
            sat='y';
        %     sat=input('are you satisfied with fit? y/n [y] ','s');
        %     if isempty(sat)
        %         sat='y';
        %     end
        %     if sat=='n';
        %         waterFreq=input('input new water frequency guess: ');
        %     end

        end


        FWHM2=abs(parsFit(2));
        FWHM2=FWHM2*(42.577*temp{ss}.Bo);  %Assumes Proton.
    else
        FWHM1 = 0;
        FWHM2 = 0;
    end

    FWHM(ss)=mean([FWHM1 FWHM2]);
end
