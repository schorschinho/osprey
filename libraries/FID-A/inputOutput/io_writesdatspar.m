%io_writesdatspar.m
%Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% RF=io_writesdatspar(in,outfile);
% 
% DESCRIPTION:
% Takes MRS jeader data in matlab structure format and writes it to a pair 
% of Philips SDAT/SPAR files that can be read by LCModel, jMRUI etc.
% 
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output SDAT/SPAR files.
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function are files in SDAT/SPAR format. 

function RF=io_writesdatspar(in,outfile);
%function RF=io_writesdatspar(in,outfile);

if ~in.flags.averaged
    disp('WARNING:  Signals must be averaged first');
end

if ~in.flags.addedrcvrs
    error('ERROR:  receiver channels must be combined first');
end

% Create dummy copy to return
RF = in;

%% Create SDAT file
% Create dummy array to hold FIDs
fids = in.fids;
if length(in.sz) == 2
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2)]);
elseif length(in.sz) == 3
    reshape_data = reshape(fids, [1 in.sz(1) in.sz(2) in.sz(3)]);
end
data=[real(reshape_data);imag(reshape_data)];

% Open a file with the designated name, and write the FID
sdatFile    = [outfile '.SDAT'];
fid         = fopen(sdatFile,'w','ieee-le');
status      = fwriteVAXG(fid,data,'float32');
fclose(fid);

%% Create SPAR file

% Open a file with the designated name, and write the header information
sparFile    = [outfile '.SPAR'];
fid         = fopen(sparFile,'w+');
fprintf(fid,'!--------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'!\n');
fprintf(fid,'\n');
fprintf(fid,'!\n');
fprintf(fid,'\n');
fprintf(fid,'!      CAUTION - Investigational device.\n');
fprintf(fid,'\n');
fprintf(fid,'!      Limited by Federal Law to investigational use.\n');
fprintf(fid,'\n');
fprintf(fid,'!\n');
fprintf(fid,'\n');
fprintf(fid,'!\n');
fprintf(fid,'\n');
fprintf(fid,'!      GYROSCAN spectro parameter file.\n');
fprintf(fid,'\n');
fprintf(fid,'!      Last revised 05-July-2007.\n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'!   This file contains time domain data in the spectral dimension.\n');
fprintf(fid,'\n');
fprintf(fid,'!   S15/ACS: set of *.SPAR and *.SDAT files is created, (dataformat: VAX CPX floats)\n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'examination_name : \n');
fprintf(fid,'\n');
fprintf(fid,'scan_id : \n');
fprintf(fid,'\n');
fprintf(fid,'scan_date : \n');
fprintf(fid,'\n');
fprintf(fid,'patient_name : \n');
fprintf(fid,'\n');
fprintf(fid,'patient_birth_date : \n');
fprintf(fid,'\n');
fprintf(fid,'patient_position : "head_first"\n');
fprintf(fid,'\n');
fprintf(fid,'patient_orientation : "supine"\n');
fprintf(fid,'\n');
fprintf(fid,'samples : %i\n', in.sz(1));
fprintf(fid,'\n');
fprintf(fid,'rows : %i\n', in.sz(2));
fprintf(fid,'\n');
fprintf(fid,'synthesizer_frequency : %i\n', in.txfrq);
fprintf(fid,'\n');
fprintf(fid,'offset_frequency : 0\n');
fprintf(fid,'\n');
fprintf(fid,'sample_frequency : %i\n', in.spectralwidth);
fprintf(fid,'\n');
fprintf(fid,'echo_nr : 1\n');
fprintf(fid,'\n');
fprintf(fid,'mix_number : 1\n');
fprintf(fid,'\n');
fprintf(fid,'nucleus : 1H\n');
fprintf(fid,'\n');
fprintf(fid,'t0_mu1_direction : 0\n');
fprintf(fid,'\n');
fprintf(fid,'echo_time : %i\n', in.te);
fprintf(fid,'\n');
fprintf(fid,'repetition_time : %i\n', in.tr);
fprintf(fid,'\n');
fprintf(fid,'averages : %i\n', in.averages);
fprintf(fid,'\n');
fprintf(fid,'volume_selection_enable : "yes"\n');
fprintf(fid,'\n');
fprintf(fid,'volumes : 1\n');
fprintf(fid,'\n');
fprintf(fid,'ap_size : %4.8f\n', in.geometry.size.ap);
fprintf(fid,'\n');
fprintf(fid,'lr_size : %4.8f\n', in.geometry.size.lr);
fprintf(fid,'\n');
fprintf(fid,'cc_size : %4.8f\n', in.geometry.size.cc);
fprintf(fid,'\n');
fprintf(fid,'ap_off_center : %4.8f\n', in.geometry.pos.ap);
fprintf(fid,'\n');
fprintf(fid,'lr_off_center : %4.8f\n', in.geometry.pos.lr);
fprintf(fid,'\n');
fprintf(fid,'cc_off_center : %4.8f\n', in.geometry.pos.cc);
fprintf(fid,'\n');
fprintf(fid,'ap_angulation : %4.8f\n', in.geometry.rot.ap);
fprintf(fid,'\n');
fprintf(fid,'lr_angulation : %4.8f\n', in.geometry.rot.lr);
fprintf(fid,'\n');
fprintf(fid,'cc_angulation : %4.8f\n', in.geometry.rot.cc);
fprintf(fid,'\n');
fprintf(fid,'volume_selection_method : 1\n');
fprintf(fid,'\n');
fprintf(fid,'phase_encoding_enable : "no"\n');
fprintf(fid,'\n');
fprintf(fid,'t1_measurement_enable : "no"\n');
fprintf(fid,'\n');
fprintf(fid,'t2_measurement_enable : "no"\n');
fprintf(fid,'\n');
fprintf(fid,'time_series_enable : "no"\n');
fprintf(fid,'\n');
fprintf(fid,'image_plane_slice_thickness : 0\n');
fprintf(fid,'\n');
fprintf(fid,'slice_distance : 0\n');
fprintf(fid,'\n');
fprintf(fid,'nr_of_slices_for_multislice : 1\n');
fprintf(fid,'\n');
fprintf(fid,'Spec.image in plane transf : "plusB-minA"\n');
fprintf(fid,'\n');
fprintf(fid,'!---------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'spec_data_type : cf\n');
fprintf(fid,'\n');
fprintf(fid,'spec_sample_extension :[V]\n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------\n');
fprintf(fid,'\n');
fprintf(fid,'! Column parameters \n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------\n');
fprintf(fid,'\n');
fprintf(fid,'spec_num_col : %i\n', in.sz(1));
fprintf(fid,'\n');
fprintf(fid,'spec_col_lower_val : %4.6f\n', -in.spectralwidth/2);
fprintf(fid,'\n');
fprintf(fid,'spec_col_upper_val : %4.6f\n', -in.spectralwidth/2 + in.sz(1)/in.spectralwidth);
fprintf(fid,'\n');
fprintf(fid,'spec_col_extension :[sec]\n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------\n');
fprintf(fid,'\n');
fprintf(fid,'! Row parameters \n');
fprintf(fid,'\n');
fprintf(fid,'!--------------------\n');
fprintf(fid,'\n');
fprintf(fid,'spec_num_row : %i\n', in.averages);
fprintf(fid,'\n');
fprintf(fid,'spec_row_lower_val : %i\n', 1);
fprintf(fid,'\n');
fprintf(fid,'spec_row_upper_val : %i\n', in.averages);
fprintf(fid,'\n');
fprintf(fid,'spec_row_extension :[index]\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'! Extra parameters in order to make data transfer \n');
fprintf(fid,'\n');
fprintf(fid,'! possible between S15/ACS and SUNSPEC: \n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'num_dimensions : 2\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'dim1_ext       :[sec]\n');
fprintf(fid,'\n');
fprintf(fid,'dim1_pnts : %i\n', in.sz(1));
fprintf(fid,'\n');
fprintf(fid,'dim1_low_val : %4.6f\n', -in.spectralwidth/2);
fprintf(fid,'\n');
fprintf(fid,'dim1_step : %4.8f\n', in.dwelltime);
fprintf(fid,'\n');
fprintf(fid,'dim1_direction : mu1\n');
fprintf(fid,'\n');
fprintf(fid,'dim1_t0_point : 0\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_ext :[num]\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_pnts : 1\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_low_val : 1.000000\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_step : 1.000000\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_direction : x\n');
fprintf(fid,'\n');
fprintf(fid,'dim2_t0_point : 50\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_ext :[num]\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_pnts : 1\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_low_val : 1.000000\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_step : 1.000000\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_direction : y\n');
fprintf(fid,'\n');
fprintf(fid,'dim3_t0_point : 50\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'!   Additional parameters\n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'echo_acquisition : ECHO\n');
fprintf(fid,'\n');
fprintf(fid,'TSI_factor : 1\n');
fprintf(fid,'\n');
fprintf(fid,'spectrum_echo_time : %i\n', in.te);
fprintf(fid,'\n');
fprintf(fid,'spectrum_inversion_time : 0\n');
fprintf(fid,'\n');
fprintf(fid,'image_chemical_shift : 0\n');
fprintf(fid,'\n');
fprintf(fid,'resp_motion_comp_technique : NONE\n');
fprintf(fid,'\n');
fprintf(fid,'de_coupling : NO\n');
fprintf(fid,'\n');
fprintf(fid,'equipment_sw_verions : 5.4.0 ; .4.0 ;\n');
fprintf(fid,'\n');
fprintf(fid,'placeholder1 : \n');
fprintf(fid,'\n');
fprintf(fid,'placeholder2 : \n');
fprintf(fid,'\n');
fprintf(fid,'!-------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fclose(fid);