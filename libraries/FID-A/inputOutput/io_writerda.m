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
if length(in.sz) == 1
    reshape_data = fids;
elseif length(in.sz) == 2
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2)]);
elseif length(in.sz) == 3
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2) in.sz(3)]);
end
data=[real(reshape_data);-imag(reshape_data)];

%% Create RDA file

% Open a file with the designated name, and write the header information
sparFile    = [outfile '.RDA'];
fid         = fopen(sparFile,'wt+');
fprintf(fid,'>>> Begin of header <<<\r\n');
fprintf(fid,'PatientName: \r\n');
fprintf(fid,'PatientID: \r\n');
fprintf(fid,'PatientSex: \r\n');
fprintf(fid,'PatientBirthDate: \r\n');
fprintf(fid,'StudyDate: \r\n');
fprintf(fid,'StudyTime: \r\n');
fprintf(fid,'StudyDescription: \r\n');
fprintf(fid,'PatientAge: \r\n');
fprintf(fid,'PatientWeight: \r\n');
fprintf(fid,'SeriesDate: \r\n');
fprintf(fid,'SeriesTime: \r\n');
fprintf(fid,'SeriesDescription: \r\n');
fprintf(fid,'ProtocolName: \r\n');
fprintf(fid,'PatientPosition: \r\n');
fprintf(fid,'SeriesNumber: \r\n');
fprintf(fid,'InstitutionName: \r\n');
fprintf(fid,'StationName: \r\n');
fprintf(fid,'ModelName: \r\n');
fprintf(fid,'DeviceSerialNumber: \r\n');
fprintf(fid,'SoftwareVersion[0]: \r\n');
fprintf(fid,'InstanceDate: \r\n');
fprintf(fid,'InstanceTime: \r\n');
fprintf(fid,'InstanceNumber: \r\n');
fprintf(fid,'InstanceComments: \r\n');
fprintf(fid,'AcquisitionNumber: \r\n');
fprintf(fid,'SequenceName: \r\n');
fprintf(fid,'SequenceDescription: \r\n');
fprintf(fid,'TR: %4.6f\r\n', in.tr);
fprintf(fid,'TE: %4.6f\r\n', in.te);
fprintf(fid,'TM: 0.000000\r\n');
fprintf(fid,'TI: 0.000000\r\n');
fprintf(fid,'DwellTime: %i\r\n', in.dwelltime*1e6);
fprintf(fid,'EchoNumber: 1\r\n');
fprintf(fid,'NumberOfAverages: %4.6f\r\n', in.averages);
fprintf(fid,'MRFrequency: %4.6f\r\n', in.txfrq*1e-6);
fprintf(fid,'MagneticFieldStrength: %4.6f\r\n', in.Bo);
fprintf(fid,'NumOfPhaseEncodingSteps: 1\r\n');
fprintf(fid,'FlipAngle: 90.000000\r\n');
fprintf(fid,'VectorSize: %i\r\n', in.sz(1));
fprintf(fid,'CSIMatrixSize[0]: 1\r\n');
fprintf(fid,'CSIMatrixSize[1]: 1\r\n');
fprintf(fid,'CSIMatrixSize[2]: 1\r\n');
fprintf(fid,'CSIMatrixSizeOfScan[0]: 1\r\n');
fprintf(fid,'CSIMatrixSizeOfScan[1]: 1\r\n');
fprintf(fid,'CSIMatrixSizeOfScan[2]: 1\r\n');
fprintf(fid,'CSIGridShift[0]: 0\r\n');
fprintf(fid,'CSIGridShift[1]: 0\r\n');
fprintf(fid,'CSIGridShift[2]: 0\r\n');
fprintf(fid,'HammingFilter: Off\r\n');
fprintf(fid,'FrequencyCorrection: NO\r\n');
fprintf(fid,'TransmitCoil: Body\r\n');
fprintf(fid,'TransmitRefAmplitude[1H]: \r\n'); % to read into extra Siemens parameter field
fprintf(fid,'SliceThickness: \r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[0]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[1]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PositionVector[2]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[0]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[1]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'RowVector[2]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[0]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[1]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'ColumnVector[2]: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'VOIPositionSag: %4.6f\r\n', in.geometry.pos.PosSag);
fprintf(fid,'VOIPositionCor: %4.6f\r\n', in.geometry.pos.PosCor);
fprintf(fid,'VOIPositionTra: %4.6f\r\n', in.geometry.pos.PosTra);
fprintf(fid,'VOIThickness: %4.6f\r\n', in.geometry.size.VoIThickness);
fprintf(fid,'VOIPhaseFOV: %4.6f\r\n', in.geometry.size.VoI_PeFOV);
fprintf(fid,'VOIReadoutFOV: %4.6f\r\n', in.geometry.size.VoI_RoFOV);
fprintf(fid,'VOINormalSag: %4.6f\r\n', in.geometry.rot.NormSag);
fprintf(fid,'VOINormalCor: %4.6f\r\n', in.geometry.rot.NormCor);
fprintf(fid,'VOIRotationInPlane: %4.6f\r\n', in.geometry.rot.VoI_InPlaneRot);
fprintf(fid,'FoVHeight: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'FoVWidth: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'FoV3D: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PercentOfRectFoV: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'NumberOfRows: 1\r\n');
fprintf(fid,'NumberOfColumns: 1\r\n');
fprintf(fid,'NumberOf3DParts: 1\r\n');
fprintf(fid,'PixelSpacingRow: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PixelSpacingCol: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'PixelSpacing3D: 0\r\n'); % to read into extra Siemens parameter field
fprintf(fid,'>>> End of header <<<\r\n');
fwrite(fid,data,'double');
fclose(fid);