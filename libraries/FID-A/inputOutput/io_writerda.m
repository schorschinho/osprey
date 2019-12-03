%io_writerda.m
%Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% RF=io_writerda(in,outfile);
% 
% DESCRIPTION:
% Takes MRS jeader data in matlab structure format and writes it to a pair 
% of Siemens RDA files that can be read by LCModel, jMRUI etc.
%
% The written RDA file header information will be entirely de-identified,
% ie fields containing PHI are not filled.
% 
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output RDA file.
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in RDA format. 

function RF=io_writerda(in,outfile);
%function RF=io_writerda(in,outfile);

if ~in.flags.averaged
    disp('WARNING:  Signals must be averaged first');
end

if ~in.flags.addedrcvrs
    error('ERROR:  receiver channels must be combined first');
end

% Create dummy copy to return
RF = in;

% Create dummy array to hold FIDs
fids = in.fids;
if length(in.sz) == 2
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2)]);
elseif length(in.sz) == 3
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2) in.sz(3)]);
end
data=[real(reshape_data);-imag(reshape_data)];

%% Create RDA file

% Open a file with the designated name, and write the header information
sparFile    = [outfile '.RDA'];
fid         = fopen(sparFile,'w+');
fprintf(fid,'>>> Begin of header <<<\n');
fprintf(fid,'PatientName: \n');
fprintf(fid,'PatientID: \n');
fprintf(fid,'PatientSex: \n');
fprintf(fid,'PatientBirthDate: \n');
fprintf(fid,'StudyDate: \n');
fprintf(fid,'StudyTime: \n');
fprintf(fid,'StudyDescription: \n');
fprintf(fid,'PatientAge: \n');
fprintf(fid,'PatientWeight: \n');
fprintf(fid,'SeriesDate: \n');
fprintf(fid,'SeriesTime: \n');
fprintf(fid,'SeriesDescription: \n');
fprintf(fid,'ProtocolName: \n');
fprintf(fid,'PatientPosition: \n');
fprintf(fid,'SeriesNumber: \n');
fprintf(fid,'InstitutionName: \n');
fprintf(fid,'StationName: \n');
fprintf(fid,'ModelName: \n');
fprintf(fid,'DeviceSerialNumber: \n');
fprintf(fid,'SoftwareVersion[0]: \n');
fprintf(fid,'InstanceDate: \n');
fprintf(fid,'InstanceTime: \n');
fprintf(fid,'InstanceNumber: \n');
fprintf(fid,'InstanceComments: \n');
fprintf(fid,'AcquisitionNumber: \n');
fprintf(fid,'SequenceName: \n');
fprintf(fid,'SequenceDescription: \n');
fprintf(fid,'TR: %4.6f\n', in.tr);
fprintf(fid,'TE: %4.6f\n', in.te);
fprintf(fid,'TM: 0.000000\n');
fprintf(fid,'TI: 0.000000\n');
fprintf(fid,'DwellTime: %i\n', in.dwelltime*1e6);
fprintf(fid,'EchoNumber: 1\n');
fprintf(fid,'NumberOfAverages: %4.6f\n', in.sz(2));
fprintf(fid,'MRFrequency: %4.6f\n', in.txfrq*1e-6);
fprintf(fid,'MagneticFieldStrength: %4.6f\n', in.Bo);
fprintf(fid,'NumOfPhaseEncodingSteps: 1\n');
fprintf(fid,'FlipAngle: 90.000000\n');
fprintf(fid,'VectorSize: %i\n', in.sz(1));
fprintf(fid,'CSIMatrixSize[0]: 1\n');
fprintf(fid,'CSIMatrixSize[1]: 1\n');
fprintf(fid,'CSIMatrixSize[2]: 1\n');
fprintf(fid,'CSIMatrixSizeOfScan[0]: 1\n');
fprintf(fid,'CSIMatrixSizeOfScan[1]: 1\n');
fprintf(fid,'CSIMatrixSizeOfScan[2]: 1\n');
fprintf(fid,'CSIGridShift[0]: 0\n');
fprintf(fid,'CSIGridShift[1]: 0\n');
fprintf(fid,'CSIGridShift[2]: 0\n');
fprintf(fid,'HammingFilter: Off\n');
fprintf(fid,'FrequencyCorrection: NO\n');
fprintf(fid,'TransmitCoil: Body\n');
fprintf(fid,'TransmitRefAmplitude[1H]: \n'); % to read into extra Siemens parameter field
fprintf(fid,'SliceThickness: \n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[0]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[1]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[2]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[0]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[1]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[2]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[0]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[1]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[2]: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'VOIPositionSag: %4.6f\n', in.geometry.pos.PosSag);
fprintf(fid,'VOIPositionCor: %4.6f\n', in.geometry.pos.PosCor);
fprintf(fid,'VOIPositionTra: %4.6f\n', in.geometry.pos.PosTra);
fprintf(fid,'VOIThickness: %4.6f\n', in.geometry.size.VoIThickness);
fprintf(fid,'VOIPhaseFOV: %4.6f\n', in.geometry.size.VoI_PeFOV);
fprintf(fid,'VOIReadoutFOV: %4.6f\n', in.geometry.size.VoI_RoFOV);
fprintf(fid,'VOINormalSag: %4.6f\n', in.geometry.rot.NormSag);
fprintf(fid,'VOINormalCor: %4.6f\n', in.geometry.rot.NormCor);
fprintf(fid,'VOIRotationInPlane: %4.6f\n', in.geometry.rot.VoI_InPlaneRot);
fprintf(fid,'FoVHeight: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'FoVWidth: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'FoV3D: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'PercentOfRectFoV: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'NumberOfRows: 1\n');
fprintf(fid,'NumberOfColumns: 1\n');
fprintf(fid,'NumberOf3DParts: 1\n');
fprintf(fid,'PixelSpacingRow: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'PixelSpacingCol: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'PixelSpacing3D: 0\n'); % to read into extra Siemens parameter field
fprintf(fid,'>>> End of header <<<\n');
fwrite(fid,data,'double');
fclose(fid);