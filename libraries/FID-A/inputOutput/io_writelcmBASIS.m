%% io_writelcmBASIS
%   This function creates a LCM usable .BASIS file from an Osprey Basisset.
%   MMs are not included in the output
%
%
%   USAGE:
%       RF=io_writelcmBASIS(in,outfile,vendor,resample);
%
%   INPUT:      in      = Osprey BASIS file
%               outfile = path and name of the LCModel .BASIS file
%               vendor  = String with Vendor name for consistent naming 
%               SEQ     = name of the sequence
%
%   OUTPUT:     RF is unused, but .BASIS file is created
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-02-11: First version of the code.

function RF=io_writelcmBASIS(in,outfile,vendor,SEQ);

metabList = fit_createMetabList({'default'});
% Create the modified basis set without macro molecules 
basisSet = fit_selectMetabs(in, metabList, 0);

Bo=basisSet.Bo;
HZPPPM=42.577*Bo;
FWHMBA = basisSet.linewidth/HZPPPM;
ECHOT = basisSet.te;
  
BADELT=basisSet.dwelltime;
NDATAB= basisSet.sz(1);

XTRASH = 0;

%write to txt file
fid=fopen(outfile,'w+');
fprintf(fid,' $SEQPAR');
fprintf(fid,'\n FWHMBA = %5.6f,',FWHMBA);
fprintf(fid,'\n HZPPPM = %5.6f,',HZPPPM);
fprintf(fid,'\n ECHOT = %2.2f,',ECHOT);
fprintf(fid,'\n SEQ = ''%s''',SEQ);
fprintf(fid,'\n $END');
fprintf(fid,'\n $BASIS1');
fprintf(fid,'\n IDBASI = ''%s'',',[vendor ' ' SEQ ' ' num2str(ECHOT) ' ms Osprey']);
fprintf(fid,'\n FMTBAS = ''(2E15.6)'',');
fprintf(fid,'\n BADELT = %5.6f,',BADELT);
fprintf(fid,'\n NDATAB = %i', NDATAB);
fprintf(fid,'\n $END\n');
for i = 1 : basisSet.nMets
    if ~strcmp(basisSet.name{i}, 'CrCH2') && ~strcmp(basisSet.name{i}, 'H2O')
        RF = shift_centerFreq(basisSet,i);
        fprintf(fid,' $NMUSED');
        fprintf(fid,'\n XTRASH = %2.2f',XTRASH);
        fprintf(fid,'\n $END');
        fprintf(fid,'\n $BASIS');
        fprintf(fid,'\n ID = ''%s'',',basisSet.name{i});
        fprintf(fid,'\n METABO = ''%s'',',basisSet.name{i});
        fprintf(fid,'\n CONC = 1.,');
        fprintf(fid,'\n TRAMP = 1.,');
        fprintf(fid,'\n VOLUME = 1.,');
        fprintf(fid,'\n ISHIFT = 0');
        fprintf(fid,'\n $END\n');
        fprintf(fid,' %7.6e  %7.6e\n',RF');
    end
end

fclose(fid);
end

function [RF] = shift_centerFreq(data_struct,idx)

    t=repmat(data_struct.t',[1 data_struct.sz(2:end,1)]);
    hzpppm = data_struct.Bo*42.577;
    f = (4.68-data_struct.centerFreq)*hzpppm;
    fids = data_struct.fids(:,idx);
    fids=fids.*exp(-1i*t*f*2*pi);
    %Take the complex conjugate becuase the sense of rotation in LCModel seems to
    %be opposite to that used in FID-A.
    fids = conj(fids);
    specs=(fft(fids,[],data_struct.dims.t));
    RF=zeros(length(specs(:)),2);
    RF(:,1)=real(specs(:));
    RF(:,2)=imag(specs(:));

end