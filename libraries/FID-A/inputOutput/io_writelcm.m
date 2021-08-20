%io_writelcm.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=io_writelcm(in,outfile,te);
% 
% DESCRIPTION:
% Takes MRS data in matlab structure format and writes it to a text file
% that can be read by LCModel.
% 
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output text file.
% te         = Echo time of acquisition (in ms).
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in LCModel raw format. 

function RF=io_writelcm(in,outfile,te);
%function RF=writelcm(in,outfile,te);

% if in.flags.isISIS
%     error('ERROR:  Must make subspecs first');
% end

if ~in.flags.averaged
    disp('WARNING:  Signals must be averaged first');
end

if ~in.flags.addedrcvrs
    error('ERROR:  receiver channels must be combined first');
end


Bo=in.Bo;
hzppm=42.577*Bo;
dwelltime=in.dwelltime;

% Calculate the voxel volume from certain fields
if isfield(in.geometry.size,'ap') % Philips
    % Philips has names attached to the fields
    vol = in.geometry.size.ap * in.geometry.size.lr * in.geometry.size.cc;
elseif isfield(in.geometry.size,'VoIThickness') % Siemens
    % Siemens has other names attached to the fields
    vol = in.geometry.size.VoI_RoFOV * in.geometry.size.VoI_PeFOV * in.geometry.size.VoIThickness;
else %GE
    % For GE data, there is currently no designation
    vol = in.geometry.size.dim1 * in.geometry.size.dim2 * in.geometry.size.dim3;
end
% Avoid zero volume
if vol == 0
    vol = 1000.0;
end


RF=zeros(in.sz(in.dims.t),2);
RF(:,1)=real(in.fids(:,1));
RF(:,2)=-imag(in.fids(:,1));


%write to txt file for jmrui
fid=fopen(outfile,'w+');
fprintf(fid,' $SEQPAR');
%fprintf(fid,'\n\nFilename: %s' ,outfile);
%fprintf(fid,'\n\nPointsInDataset: %i',length(data_struct.fids));
%fprintf(fid,'\nDatasetsInFile: %i',datsets);
%fprintf(fid,'\nSamplingInterval: %1.0E',data_struct.dwelltime*1000);
%fprintf(fid,'\nZeroOrderPhase: %1.0E',zop);
%fprintf(fid,'\nBeginTime: %1.0E',t0);
%fprintf(fid,'\nTransmitterFrequency: %2.4E',data_struct.txfrq);
%fprintf(fid,'\nMagneticField: %2.1E',Bo);
%fprintf(fid,'\nTypeOfNucleus: %1.0E',Nuc);
%fprintf(fid,'\nNameOfPatient: %s',PatName);
%fprintf(fid,'\nDateOfExperiment: %i',data_struct.date);
%fprintf(fid,'\nSpectrometer: %s',scanner);
%fprintf(fid,'\nAdditionalInfo: %s\n\n\n',addinfo);
%fprintf(fid,'Signal and FFT\n');
%fprintf(fid,'sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n');
%fprintf(fid,'Signal 1 out of %i in file\n',datsets);
fprintf(fid,'\n ECHOT= %2.2f',te);
if isfield(in, 'seq')
    fprintf(fid,'\n SEQ= ''%s''', in.seq);
end
fprintf(fid,'\n HZPPPM= %5.6f',in.txfrq/1e6);
fprintf(fid,'\n NUNFIL= %i',in.sz(1));
fprintf(fid,'\n DELTAT= %5.6f' ,dwelltime);
fprintf(fid,'\n $END');
fprintf(fid,'\n $NMID');
fprintf(fid,'\n id=''ANONYMOUS '', fmtdat=''(2E15.6)''');
fprintf(fid,'\n volume=%2.1f', vol/1e3);
fprintf(fid,'\n tramp=1.0');
fprintf(fid,'\n $END\n');
fprintf(fid,'  % 7.6e  % 7.6e\n',RF');
fclose(fid);
