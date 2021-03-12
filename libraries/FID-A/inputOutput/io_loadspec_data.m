%[out,out_w]=io_loadspec_data(filename,subspecs,kk,statFile);
%
% USAGE:
% [out,out_w]=io_loadspec_data(ffilename,subspecs,kk,statFile);
% 
% DESCRIPTION:
% Reads in philips MRS data (.data and .list files) using code adapted from 
% PhilipsRead_data.m, provided as part of the Gannet software package by 
% Richard Edden (gabamrs.blogspot.com).
% 
% op_loadspec_data outputs the data in structure format, with fields 
% corresponding to time scale, fids, frequency scale, spectra, and header 
% fields containing information about the acquisition.  The resulting 
% matlab structure can be operated on by the other functions in this MRS 
% toolbox.  NOTE:  Since the Gannet code is geared towards edited GABA MRS 
% data, this code may not be general enough to handle all types of MRS data.  
% Suggestions are most welcome.
% 
% INPUTS:
% filename   = filename of Philips .data file to be loaded.
% subspecs   = number of subspectra in the data (from spectral editing, ISIS, etc.)
% kk         = file number in container (only used if sequence parameters are stored in statFile.csv file).
% statFile   = Path to statFile.csv file(Optional if no SPAR file is available).
%
% OUTPUTS:
% out        = Input water suppressed dataset in FID-A structure format.
% out_w      = Input water reference dataset in FID-A structure format. 
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-10-03)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2020-10-02: First version of the code.

function [out,out_w]=io_loadspec_data(filename,subspecs,kk,statFile);
%% Parse Input
if nargin < 4
    statFile = [];
end

%% Extract information from SPAR files that is not in DATA/LIST
hasSPAR = 0;
[~,~,ext] = fileparts(filename);
% Get the appropriate raw_act.spar
spar_file = strrep(filename,ext, '_raw_act.spar'); % it's automatically case-insensitive

if ~isfile(spar_file)
    spar_file = strrep(filename,ext, '_raw_act.SPAR'); % it's automatically case-insensitive
end
if ~isfile(spar_file)    
    if ~isempty(statFile) % Has stat csv file
        statCSV = readtable(statFile, 'Delimiter', ',','ReadVariableNames',1); % Load it
        name = statCSV.Properties.VariableNames;
        tr_idx = find(strcmp(name,'tr'));
        if isempty(tr_idx)
            tr_idx = find(strcmp(name,'TR'));
        end
        te_idx = find(strcmp(name,'te'));        
        if isempty(te_idx)
            te_idx = find(strcmp(name,'TE'));
        end
        sw_idx = find(strcmp(name,'sw'));        
        if isempty(sw_idx)
            sw_idx = find(strcmp(name,'SW'));
        end 
        Larmor_idx = find(strcmp(name,'Larmor')); 
        date_idx = find(strcmp(name,'date'));    
        seq_idx = find(strcmp(name,'seq')); 
        if ~isempty(tr_idx) && ~isempty(te_idx) && ~isempty(sw_idx) && ~isempty(Larmor_idx)
            tr = statCSV{kk,tr_idx};
            te = statCSV{kk,te_idx}; 
            sw = statCSV{kk,sw_idx}; 
            Larmor = statCSV{kk,Larmor_idx}; 
            if ~isempty(date_idx)
                date = statCSV{kk,date_idx}{1};
            else
                date = '';
            end
            if ~isempty(seq_idx)
                seq = statCSV{kk,seq_idx}{1};
            else
                seq = 'PRESS'; % Let's assume it is a PRESS sequence 
            end
        else
            msg = 'You need a SPAR file that has the same name as the .data/.list file with the extension _raw_act.spar. Otherwise you can supply TR, TE, SW (spectralwidth), Larmor (larmor frequency in MHz), and seq (localization) in the stat.csv file.';
            error(msg);    
        end
    else
        msg = 'You need a SPAR file that has the same name as the .data/.list file with the extension _raw_act.spar. Otherwise you can supply TR, TE, SW (spectralwidth), and Larmor (larmor frequency in MHz), and seq (localization) in the stat.csv file.';
        error(msg);          
    end
else
    hasSPAR = 1;    
end

if hasSPAR
    % Open spar file and read parameters
    sparname = fopen(spar_file,'r');
    sparheader = textscan(sparname, '%s');
    sparheader = sparheader{1};
    sparidx=find(ismember(sparheader, 'repetition_time')==1); % TR
    tr = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'echo_time')==1); % TE
    te = str2double(sparheader{sparidx+2}); % 
    sparidx=find(ismember(sparheader, 'synthesizer_frequency')==1); % F0
    Larmor = str2double(sparheader{sparidx+2})/1e6;
    sparidx=find(ismember(sparheader, 'sample_frequency')==1); % readout bandwidth
    sw = str2double(sparheader{sparidx+2});

    % Read voxel geometry information.
    % THIS IS IN THE ORDER LR-AP-FH
    sparidx=find(ismember(sparheader, 'scan_date')==1); % voxel size 
    date = sparheader{sparidx+2};
    sparidx=find(ismember(sparheader, 'scan_id')==1); % voxel size 
    seq = sparheader{sparidx+2};
    sparidx=find(ismember(sparheader, 'lr_size')==1); % voxel size 
    geometry.size.lr = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'ap_size')==1);
    geometry.size.ap = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'cc_size')==1);
    geometry.size.cc = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'lr_off_center')==1); % voxel center offset
    geometry.pos.lr = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'ap_off_center')==1);
    geometry.pos.ap = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'cc_off_center')==1);
    geometry.pos.cc = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'lr_angulation')==1); % voxel angulation (radians)
    geometry.rot.lr = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'ap_angulation')==1);
    geometry.rot.ap = str2double(sparheader{sparidx+2});
    sparidx=find(ismember(sparheader, 'cc_angulation')==1);
    geometry.rot.cc = str2double(sparheader{sparidx+2});
    fclose(sparname);
end

%read in the data using the philipsDataLoad.m (adapted from PhilipsRead_data.m)
[data] = loadRawKspace(filename);
%[FullData,WaterData]=philipsDataLoad(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine dimensions of the acquisition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine number of channels
    n_coils = length(unique(data.chan));
    
    % Determine number of mixes
    n_mixes = length(unique(data.mix));
    
    % Determine number of averages per mix
    n_averages = data.kspace_properties.number_of_signal_averages;
    n_av_edit = n_averages(1);
    if n_mixes == 2
        n_av_water = n_averages(2); % if second mix exists, it is water
    end
    % This will be the NSA as specified on the exam card for the water-suppressed mix (mix = 0)
    % and the NSA as specified on the exam card for the water-suppressed mix (mix = 1);
    
    % Determine number of dynamics per NSA. It is not stored in the dynamics
    % field, but rather in extra attributes. Could be different for different
    % software versions, need to check!
    n_dyns = data.kspace_properties.number_of_extra_attribute_1_values;
    n_dyns_edit = n_dyns(1);
    %if n_mixes == 2
    %    n_dyns_water = n_dyns(2); % if second mix exists, it is water
    %end
    % Since dynamics are set globally on the exam card, this will be the same
    % for both - it will only be the true value for the water-suppressed mix
    
    % Determine number of data points per scan
    n_points = data.kspace_properties.F_resolution(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start splitting the list of total scans into its parts:
    % Noise scans, water-suppressed scans, and water-unsuppressed scans.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Noise scans have the data type 'NOI'
    isnoise = strcmp(data.typ,'NOI');
    fids_noise = cell2mat(data.complexdata(isnoise));
    n_av_noise = size(fids_noise,1) ./ n_coils;
    
    % Separate the water-suppressed from the water-unsuppressed scans.
    % Water-suppressed scans have the data type 'STD' and mix index 0:
    isedit = strcmp(data.typ,'STD') & (data.mix == 0);
    fids_edit = cell2mat(data.complexdata(isedit));
    % Reshape the water-suppressed data
    fids_edit = reshape(fids_edit,[n_coils n_av_edit*n_dyns_edit n_points]);
    
    %Permute order of the data to match MRSCont structure
    fids_edit = permute(fids_edit,[3 1 2]);
    
    % Phase by multiplying with normalized complex conjugate of first point
%     if ~MRS_struct.p.PRIAM
%         % Don't do this when using SENSE reconstruction!
%         conj_norm = conj(fids_edit(:,:,1)) ./ abs(conj(fids_edit(:,:,1)));
%         fids_edit_ph = fids_edit .* repmat(conj_norm, [1 1 n_points]);
%     else
%         fids_edit_ph = fids_edit;
%     end
    fids_edit_ph = fids_edit;

    
    % Water-unsuppressed scans have the data type 'STD' and mix index 1:
    if n_mixes == 2
        iswater = strcmp(data.typ,'STD') & (data.mix == 1);
        fids_water = cell2mat(data.complexdata(iswater));
        
        % Reshape the water-unsuppressed data
        fids_water = reshape(fids_water,[n_coils n_av_water n_points]);
        
        % Phase by multiplying with normalized complex conjugate of first point
%         if ~MRS_struct.p.PRIAM
%             % Don't do this when using SENSE reconstruction!
%             conj_norm = conj(fids_water(:,:,1)) ./ abs(conj(fids_water(:,:,1)));
%             fids_water_ph = fids_water .* repmat(conj_norm, [1 1 n_points]);
%         else
%             fids_water_ph = fids_water;
%         end
        %Permute order of the data to match MRSCont structure
        fids_water = permute(fids_water,[3 1 2]);
        fids_water_ph = fids_water;        


    end
    
    % Perform coil combination if not PRIAM
%     if ~MRS_struct.p.PRIAM
%         if n_mixes == 2
%             firstpoint_water = fids_water_ph(:,:,1);
%             channels_scale = squeeze(sqrt(sum(firstpoint_water .* conj(firstpoint_water),1)));
%             channels_scale = repmat(channels_scale, [size(fids_water_ph,1) 1 size(fids_water_ph,3)]);
%             fids_water_ph = fids_water_ph ./ channels_scale;
%             fids_water_ph = conj(squeeze(sum(fids_water_ph,1))).';
%             
%             MRS_struct.fids.data_water = fids_water_ph;
%             
%             channels_scale = mean(channels_scale,2);
%             channels_scale = repmat(channels_scale, [1 size(fids_edit_ph,2) 1]);
%             fids_edit_ph = fids_edit_ph ./ channels_scale;
%             fids_edit_ph = conj(squeeze(sum(fids_edit_ph,1))).';
%             
%             MRS_struct.fids.data = fids_edit_ph;
%         end
%     end
    

%As far as I can tell, the data that comes out of the philipsDataLoad
%function is normally a N x Navgs matrix.  The coils have already been 
%combined. The Navgs dimension contains all the subspectra, so we will 
%split them now.  Note, that in the data-list format that I have seen, the
%edit-OFF subspectra appear in elements [1 2 5 6 9 10 13 14...] and the
%edit-ON subspectra appear in the elements [3 4 7 8 11 12 15 16...].  Other sequences
%may result in a different subspecs order, but for now we will separate the 
%subspectra in this way.
%If the data has multiple subspectra 
if subspecs == 2
    %First make an vector that holds the indices of the ON subspectra:
    totalAvgs=size(fids_edit_ph,3);
    OFFindices=[1:2:totalAvgs]-mod([0:(totalAvgs/2)-1],2);
    ONindices=[2:2:totalAvgs]+mod([1:totalAvgs/2],2);
    %Now split the subspectra out of the "averages" dimension:
    dummy(:,:,:,1)=fids_edit_ph(:,:,OFFindices);
    dummy(:,:,:,2)=fids_edit_ph(:,:,ONindices);
    fids_edit_ph = dummy;
else if subspecs == 4     
    %First make an vector that holds the indices of the multiplexed subspectra:
        totalAvgs=size(fids_edit_ph,3);
        Aindices=[1:4:totalAvgs];
        Bindices=[2:4:totalAvgs];
        Cindices=[3:4:totalAvgs];
        Dindices=[4:4:totalAvgs];
    %Now split the subspectra out of the "averages" dimension:
        dummy(:,:,:,1)=fids_edit_ph(:,:,Aindices);
        dummy(:,:,:,2)=fids_edit_ph(:,:,Bindices);
        dummy(:,:,:,3)=fids_edit_ph(:,:,Cindices);
        dummy(:,:,:,4)=fids_edit_ph(:,:,Dindices);
        fids_edit_ph = dummy;
    end
end

fids=squeeze(fids_edit_ph);
if n_mixes == 2
    fids_w=squeeze(fids_water_ph);
end

sz=size(fids);
if n_mixes == 2
    sz_w=size(fids_w);
end

%Find the magnetic field strength:
Bo=Larmor/42.577;

%Find the number of averages:
% Naverages=size(fids,3)*size(fids,4);
% Naverages_w=size(fids_w,3)*size(fids_w,4);

%In Philips data/list format, coil channels have already been combined:
% Ncoils=1;
% Ncoils_w=1;

%Now create a record of the dimensions of the data array.  
dims.t=1;
dims.coils=2;
dims.averages=3;
if subspecs>1
    dims.subSpecs=4;
else
    dims.subSpecs=0;
end
dims.extras=0;

if n_mixes == 2
    dims_w.t=1;
    dims_w.coils=2;
    dims_w.averages=3;
    dims_w.subSpecs=0;
    dims_w.extras=0;
end

specs=fftshift(ifft(fids,[],dims.t),dims.t);
if n_mixes == 2
    specs_w=fftshift(ifft(fids_w,[],dims_w.t),dims_w.t);
end


%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
spectralwidth=sw;
dwelltime=1/spectralwidth;
    
%Get TxFrq
txfrq=Larmor*1e6;


%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
%FOR WATER SUPPRESSED DATA:
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=1;
    end
end

%FOR WATER UNSUPPRESSED DATA:
if n_mixes == 2
    if dims_w.subSpecs ~=0
        if dims_w.averages~=0
            averages_w=sz_w(dims_w.averages)*sz(dims_w.subSpecs);
            rawAverages_w=averages_w;
        else
            averages_w=sz_w(dims_w.subSpecs);
            rawAverages_w=1;
        end
    else
        if dims_w.averages~=0
            averages_w=sz_w(dims_w.averages);
            rawAverages_w=averages_w;
        else
            averages_w=1;
            rawAverages_w=1;
        end
    end
end


%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
%FOR WATER SUPPRESSED DATA:
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%FOR WATER UNSUPPRESSED DATA:
if n_mixes == 2
    if dims_w.subSpecs ~=0
        subspecs_w=sz_w(dims.subSpecs);
        rawSubspecs_w=subspecs_w;
    else
        subspecs_w=1;
        rawSubspecs_w=subspecs_w;
    end
end

%****************************************************************

%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(Bo*42.577);

% Philips data assumes the center frequency to be 4.68 ppm:
centerFreq = 4.68;
ppm=ppm + centerFreq;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FOR WATER SUPPRESSED DATA
%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;
out.centerFreq = centerFreq;
if hasSPAR
    out.geometry = geometry;
end

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end


%FOR WATER UNSUPPRESSED DATA
%FILLING IN DATA STRUCTURE
if n_mixes == 2
    out_w.fids=fids_w;
    out_w.specs=specs_w;
    out_w.sz=sz_w;
    out_w.ppm=ppm;  
    out_w.t=t;    
    out_w.spectralwidth=spectralwidth;
    out_w.dwelltime=dwelltime;
    out_w.txfrq=txfrq;
    out_w.date=date;
    out_w.dims=dims_w;
    out_w.Bo=Bo;
    out_w.averages=averages_w;
    out_w.rawAverages=rawAverages_w;
    out_w.subspecs=subspecs_w;
    out_w.rawSubspecs=rawSubspecs_w;
    out_w.seq='';
    out_w.te=te;
    out_w.tr=tr;
    out_w.pointsToLeftshift=0;
    out_w.centerFreq = centerFreq;
    if hasSPAR
        out_w.geometry = geometry;
    end


    %FILLING IN THE FLAGS
    out_w.flags.writtentostruct=1;
    out_w.flags.gotparams=1;
    out_w.flags.leftshifted=0;
    out_w.flags.filtered=0;
    out_w.flags.zeropadded=0;
    out_w.flags.freqcorrected=0;
    out_w.flags.phasecorrected=0;
    out_w.flags.averaged=0;
    out_w.flags.addedrcvrs=0;
    out_w.flags.subtracted=0;
    out_w.flags.writtentotext=0;
    out_w.flags.downsampled=0;
    if out_w.dims.subSpecs==0
        out_w.flags.isISIS=0;
    else
        out_w.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
    end
end

% No water in list file
if n_mixes ~= 2
    out_w = [];
end


%DONE
