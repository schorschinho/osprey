% op_autophase.m
% Helge Zollner, The Johns Hopkins University University 2021.
% 
% MRS_FIRSTPHC applies automatic first-order phase correction of a spectrum.     
% The optimal zero-order and first-order phase corrections for a NMR spectrum  
% are determined by minimizing entropy. The objective function is constructed 
% using a Shannon-type information entropy measure.
%
% REFERENCE: 
% L. Chen, Z.Q. Weng, L.Y. Goh and M. Garland, Journal of Magnetic Resonance. 
% 158, 164–168 (2002).
%
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
%
% Copyright (c) 2013, University of Nottingham. All rights reserved. 
% USAGE:
% [out,ph]=op_autophase_ph0ph1(in,ppmmin,ppmmax,ph,dimNum);
% 
% DESCRIPTION:
% Search for the peak located between ppmmin and ppmmax, and then phase the
% spectrum so that that peak reaches the desired phase.
% 
% INPUTS:
% in         = input data in matlab structure format.
% ppmmin     = minimum of ppm search range.
% ppmmax     = maximum of ppm search range.
% ph         = desired phase value in degrees [optional.  Default=0].
% dimNum     = which subSpec dimension to use for phasing? [Only for use in
%              data with multiple subSpectra].  
%
% OUTPUTS:
% out        = Output following automatic phasing.
% phaseShift = The phase shift (in degrees) that was applied.

function [out,phShft,phc1]=op_autophase_ph0ph1(in,ppmmin,ppmmax,ph,dimNum);


if in.dims.coils>0
    error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
end
if in.dims.averages>0
    error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
end
if in.dims.extras>0
    error('ERROR:  Can not operate on data with extras dimension!  ABORTING!!');
end
if in.dims.subSpecs>0
    if nargin<5
        plot(in.ppm,in.specs);
        legend('subspec 1','subspec 2');
        dimNum=input('Input which subspectrum to use for phasing: ');
        if nargin<4
            ph=0;
        end
    end
else
    if nargin<4
        ph=0;
    end
    dimNum=1;
end

%Zeropad the data if it hasn't already been done
if ~in.flags.zeropadded
    in_zp=op_zeropad(in,10);
else
    in_zp=in;
end

%Narrow the frequency range:
in_zp=op_freqrange(in_zp,ppmmin,ppmmax);

options=optimset('TolX',1e-8,'MaxFunEvals',1e8, 'MaxIter',1e8);
phc=fminsearch(@(x) mrs_entropy(x, in_zp), [0 0], options);

phc0 = phc(1);
phc1 = phc(2);

% Let's try this first with the internal functions
% n=length(spectrum);
% ph=(phc0+phc1.*(1:n)/n); % linear phase
% spectrum_ph=mrs_rephase(spectrum, ph);

%Now phase shift the dataset so that the desired peak has the correct phase:
phShft = phc0 + ph;
out=op_addphase(in,phShft,phc1);
end

function f = mrs_entropy(x,in)
% Entropy is defined as the normalized derivative of the NMR spectral data
% ARGS :
% x = phc0 and phc1
% spectrum = a spectrum before automatic first-order phase correction 
% RETURNS : 
% f = entropy value (Using the first derivative)

    %initial parameters
    stepsize=1; 
    func_type=1;    
    

    %dephase
    phc0=x(1);
    phc1=x(2);
    
    
    phst=op_addphase(in, phc0,phc1);
    spectrum = phst.specs;
    L=length(spectrum);

    % Calculation of first derivatives 
    if (func_type == 1)
        ds1 = abs((spectrum(3:L)-spectrum(1:L-2))/(stepsize*2));
    else
        ds1 = ((spectrum(3:L)-spectrum(1:L-2))/(stepsize*2)).^2;
    end  
    p1 = ds1./sum(ds1);

    %Calculation of Entropy
    [M,K]=size(p1);
    for i=1:M
        for j=1:K
            if (p1(i,j)==0)%in case of ln(0)
               p1(i,j)=1; 
            end
        end
    end
    h1  = -p1.*log(p1);
    H1  = sum(h1);
    %Calculation of penalty
    Pfun	= 0.0;
    as      = spectrum - abs(spectrum);
    sumas   = sum(sum(as));
    if (sumas < 0)
       Pfun = Pfun + sum(sum((as./2).^2));
    end
    P = 1000*Pfun;

    % The value of objective function
    f = H1+P;

end
