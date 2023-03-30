% io_readlcmraw_basis.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_readlcmraw_basis(filename);
%
% DESCRIPTION:
% Reads entire LCModel .basis file into multiple FID-A data structures in MATLAB.
%
% INPUTS:
% filename   = filename of LCModel .basis file.
% conjugate       = apply complex conjugate
%
% OUTPUTS:
% out        = Input basis set saved as a structure in which each field is
%               an individual metabolite basis spectrum in FID-A structure 
%               format.

function [out]=io_readlcmraw_basis(filename,conjugate)

if nargin < 2
    conjugate = 1;
end

% Begin to read the data.
fid=fopen(filename);
line=fgets(fid);

%look for FWHMBA
fwhmba_index=contains(line,'FWHMBA');
while ~fwhmba_index;
    line=fgets(fid);
    fwhmba_index=contains(line,'FWHMBA');
    line = GetNumFromString(line);
end
linewidth=str2num(line);

%look for HZPPM
hzpppm_index=contains(line,'HZPPPM');
while ~hzpppm_index;
    line=fgets(fid);
    hzpppm_index=contains(line,'HZPPPM');
    line = GetNumFromString(line);
end
hzpppm=str2num(line);
Bo=hzpppm/42.577;
linewidth=linewidth*hzpppm;

%look for TE
te_index=contains(line,'ECHOT');
while ~te_index;
    line=fgets(fid);
    te_index=contains(line,'ECHOT');
    line = GetNumFromString(line);
end
te=str2num(line);

%look for spectral width
badelt_index=contains(line,'BADELT');
while ~badelt_index;
    line=fgets(fid);
    badelt_index=contains(line,'BADELT');
    line = GetNumFromString(line);
end
dwelltime=str2num(line);
spectralwidth=1/dwelltime;
fileEnd=false;

while ~feof(fid)
     %look for a center frequency
    ppmsep_index=contains(line,'PPMSEP');    
    while ~ppmsep_index;
        line=fgets(fid);
        if ischar(line)
            ppmsep_index=contains(line,'PPMSEP');
            if contains(line,'METABO')
                break
            end
        end       
    end
    if ppmsep_index
        line = GetNumFromString(line);
        centerFreq=str2num(line);
    else
        centerFreq = [];
    end
    
    %Look for the metabolite name
    metab_index=contains(line,'METABO');
    while ~metab_index;
        line=fgets(fid);
        if ischar(line)
            metab_index=contains(line,'METABO');
            if metab_index
                break
            end
        end
    end
    if strcmp(line(end-3), '''') 
        metabName=line(10:end-4);
    else
        metabName=line(10:end-3);
    end

    
    hdrEnd_index=contains(line,'$END');
    while ~hdrEnd_index;
        line=fgets(fid);
        if ischar(line)
            hdrEnd_index=contains(line,'$END');
            if hdrEnd_index
                break
            end
        end
    end
    
    line=fgets(fid);
    
    nmused_index=contains(line,'$NMUSED');
    basis_index=contains(line,'$BASIS');
    linenum=1;
    RF=[];
    % If the line is empty skip it
    while ~isempty(line) && ~nmused_index && ~basis_index && ~fileEnd
        %dataline=line(1:semicol_index-2);
        [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
        % If read failed, output the error
        if ~isempty(errmsg);
            fclose(fid);
            error('READLCMRAW_BASIS failed with read error: %s', errmsg);
        end
        % Store the read values into rf array
        RF = [ RF ; A ];
        if feof(fid)
            fileEnd=true;
        end
        linenum = linenum + 1;
        line=fgets(fid);
        if ischar(line)
            nmused_index=contains(line,'$NMUSED');
            basis_index=contains(line,'$BASIS');
        end
    end
    specs=RF(1:2:end) + 1i*RF(2:2:end);

    % GO 2022/01/24
    % LCModel uses negative BADELT values to encrypt basis sets
    % (LCModel.f source code lines 3666 and following)
    if dwelltime < 0
        dix = 1499;
        for rr = 1:length(specs)
            [randomresult, dix] = lcmodelrng(dix);
            specs(rr) = -specs(rr) .* exp(-20*randomresult + 10);
        end
    end

    if conjugate
        specs=flipud(fftshift(conj(specs),1));
    else
        specs=(fftshift(specs,1));
    end
    vectorsize=length(specs);
    sz=[vectorsize 1];
    if mod(vectorsize,2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids=ifft(ifftshift(specs,1),[],1);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids=ifft(circshift(ifftshift(specs,1),1),[],1);
    end
    f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
    ppm=f/(Bo*42.577);
    ppm=ppm+4.68;
    t=[dwelltime:dwelltime:vectorsize*dwelltime];
    txfrq=hzpppm*1e6;
    metabName = RemoveWhiteSpaces(metabName);
    if strcmp(metabName,'-CrCH2')
        metabName = 'CrCH2';
    end
    if strcmp(metabName,'2HG')
        metabName = 'bHG';
    end
    eval(['out.' metabName '.fids=fids;']);
    eval(['out.' metabName '.specs=specs;']);
    eval(['out.' metabName '.sz=[vectorsize 1 1 1];']);
    eval(['out.' metabName '.n=vectorsize;']);
    eval(['out.' metabName '.spectralwidth=abs(spectralwidth);']);
    eval(['out.' metabName '.sz=sz;']);
    eval(['out.' metabName '.Bo=Bo;']);
    eval(['out.' metabName '.te=te;']);
    eval(['out.' metabName '.tr=[];']);
    eval(['out.' metabName '.dwelltime=abs(1/spectralwidth);']);
    eval(['out.' metabName '.linewidth=linewidth;']);
    eval(['out.' metabName '.ppm=ppm;']);
    eval(['out.' metabName '.t=t;']);
    eval(['out.' metabName '.txfrq=txfrq;']);
    if ~isempty(centerFreq)
            eval(['out.' metabName '.centerFreq=centerFreq;']);
    end
    eval(['out.' metabName '.date=date;']);
    eval(['out.' metabName '.seq='''';']);
    eval(['out.' metabName '.sim='''';']);
    eval(['out.' metabName '.dims.t=1;']);
    eval(['out.' metabName '.dims.coils=0;']);
    eval(['out.' metabName '.dims.averages=0;']);
    eval(['out.' metabName '.dims.subSpecs=0;']);
    eval(['out.' metabName '.dims.extras=0;']);
    eval(['out.' metabName '.averages=1;']);
    eval(['out.' metabName '.flags.writtentostruct=1;']);
    eval(['out.' metabName '.flags.gotparams=1;']);
    eval(['out.' metabName '.flags.leftshifted=1;']);
    eval(['out.' metabName '.flags.filtered=0;']);
    eval(['out.' metabName '.flags.zeropadded=0;']);
    eval(['out.' metabName '.flags.freqcorrected=0;']);
    eval(['out.' metabName '.flags.phasecorrected=0;']);
    eval(['out.' metabName '.flags.averaged=1;']);
    eval(['out.' metabName '.flags.addedrcvrs=1;']);
    eval(['out.' metabName '.flags.subtracted=1;']);
    eval(['out.' metabName '.flags.writtentotext=1;']);
    eval(['out.' metabName '.flags.downsampled=0;']);
    eval(['out.' metabName '.flags.isFourSteps=0;']);
    
end

fclose(fid);
end

% RF=RF';
% rf(:,1)=RF(:,2)*180/pi;
% rf(:,2)=RF(:,1);
% rf(:,3)=ones(length(RF(:,1)),1);

function str3 = GetNumFromString(str)
str1 = regexprep(str,'[,;=]', ' ');
str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
end

function str = RemoveWhiteSpaces(str)
    pattern = '[ \t\n]'; % Match zero or more spaces, tabs, or newlines, followed by a double quote
    replacement = ''; % Replace the matched string with just a double quote
    str = regexprep(str, pattern, replacement);
end


% GO 2022/01/24
% LCModel uses negative BADELT values to encrypt basis sets
% (LCModel.f source code lines 3666 and following)
%++++++++++++++++ doubleprecision VERSION 2DP (MAR 1984) ++++++++++++++    4222
%  FUNCTION RANDOM.  PRODUCES A PSEUDORANDOM REAL ON THE OPEN INTERVAL      4223
%      (0.,1.).                                                             4224
%  DIX (IN doubleprecision) MUST BE INITIALIZED TO A WHOLE NUMBER          4225
%      BETWEEN 1.0D0 AND 2147483646.0D0 BEFORE THE FIRST CALL TO RANDOM       4226
%      AND NOT CHANGED BETWEEN SUCCESSIVE CALLS TO RANDOM.                  4227
%  BASED ON L. SCHRAGE, ACM TRANS. ON MATH. SOFTWARE 5, 132 (1979).         4228
%-----------------------------------------------------------------------    4229
function [randomresult,dix]=lcmodelrng(dix);

a   = 16807.0d0;
b15 = 32768.d0;
b16 = 65536.0d0;
p   = 2147483647.0d0;

 %                                                                           4231
 %  PORTABLE RANDOM NUMBER GENERATOR                                         4232
 %   USING THE RECURSION                                                     4233
 %    DIX = DIX*A MOD P                                                      4234
 %                                                                           4235
 %                                                                           4237
 %  7**5, 2**15, 2**16, 2**31-1                                              4238
                                                             
 %  GET 15 HI ORDER BITS OF DIX                                              4241
 xhi = dix./b16;
 xhi = xhi - rem(xhi,1.0d0);
 %  GET 16 LO BITS IF DIX AND FORM LO PRODUCT                                4244
 xalo =(dix-xhi.*b16).*a;
 %  GET 15 HI ORDER BITS OF LO PRODUCT                                       4246
 leftlo = xalo./b16;
 leftlo = leftlo - rem(leftlo,1.0d0);
 %  FORM THE 31 HIGHEST BITS OF FULL PRODUCT                                 4249
 fhi = xhi.*a + leftlo;
 %  GET OVERFLO PAST 31ST BIT OF FULL PRODUCT                                4251
 k = fix(fhi./b15);
 k = fix(k - rem(k,1.0d0));
 %  ASSEMBLE ALL THE PARTS AND PRESUBTRACT P                                 4254
 %   THE PARENTHESES ARE ESSENTIAL                                           4255
 dix =(((xalo-leftlo.*b16)-p)+(fhi-k.*b15).*b16) + k;
 %  ADD P BACK IN IF NECESSARY                                               4257
 if(dix < 0.0d0)
    dix = dix + p;
 end
 %  MULTIPLY BY 1/(2**31-1)                                                  4259
 randomresult = dix.*4.656612875d-10;
end %function random