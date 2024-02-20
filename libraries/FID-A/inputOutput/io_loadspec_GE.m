
%io_loadspec_GE.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [out,out_w]=io_loadspec_GE(filename,subspecs);
% 
% DESCRIPTION:
% Reads in GE P-files (.7 file) using code adapted from GERead.m, provided
% as part of the Gannet software package by Richard Edden (gabamrs.com).
% 
% op_loadspec_GE outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.  NOTE:  Since the 
% Gannet code is geared towards edited GABA MRS data, this code may not be 
% general enough to handle all types of MRS data.  Suggestions are most welcome.
% 
% INPUTS:
% filename   = filename of GE P file to be loaded.
% subspecs   = number of subspectra in the data (from spectral editing, ISIS, etc.)
%
% OUTPUTS:
% out        = output water suppressed dataset in FID-A structure format.
% out_ref    = output eddy current water reference  dataset in FID-A structure format.
% out_w      = output quantification water reference dataset in FID-A structure format. 

function [out,out_ref,out_w]=io_loadspec_GE(filename,subspecs)

%read in the data using the GELoad.m (adapted from GERead.m)
[GEout,GEout_ref,GEhdr]=GELoad(filename);

%As far as I can tell, the data that comes out of the GELoad
%function is normally a N x Navgs x Ncoils matrix.  The Navgs dimension
%contains all the subspectra, so we will split them now:
%If the data has multiple subspectra 
if subspecs == 4  %HERMES/HERCULES
    %Split the subspectra out of the "averages" dimension:
    data(:,:,:,1)=GEout(:,1:4:end,:);
    data(:,:,:,2)=GEout(:,2:4:end,:);
    data(:,:,:,3)=GEout(:,3:4:end,:);
    data(:,:,:,4)=GEout(:,4:4:end,:);
else if  subspecs==2 %MEGA   
        %Split the subspectra out of the "averages" dimension:
        data(:,:,:,1)=GEout(:,1:2:end,:);
        data(:,:,:,2)=GEout(:,2:2:end,:);
    else
        data=GEout;
    end
end

fids=squeeze(data);
fids_ref=squeeze(GEout_ref);

%swap the averages and the coils dimensions:
fids=permute(fids,[1,3,2,4]);
fids_ref=permute(fids_ref,[1,3,2]);

twoRefs=false;  %Flag to identify if there are automatic water reference scans acquired.  There are none by default.

if GEhdr.cv24 >= 16384 && subspecs == 1
    twoRefs = true;
else
end

sz=size(fids);
sz_w=size(fids_ref);

%Find the magnetic field strength:
Bo=GEhdr.Larmor/42.577;

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

dims_ref.t=1;
dims_ref.coils=2;
dims_ref.averages=3;
dims_ref.subSpecs=0;
dims_ref.extras=0;

specs=fftshift(fft(fids,[],dims.t),dims.t);
specs_w=fftshift(fft(fids_ref,[],dims_ref.t),dims_ref.t);


%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
spectralwidth=GEhdr.sw;
dwelltime=1/spectralwidth;
    
%Get TxFrq
txfrq=GEhdr.Larmor*1e6;

%Leave date blank
date='';

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
%FOR WATER SUPPRESSED DATA:
if dims.subSpecs~=0
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
if dims_ref.subSpecs~=0
    if dims_ref.averages~=0
        averages_w=sz_w(dims_ref.averages)*sz_w(dims_ref.subSpecs);
        rawAverages_w=averages_w;
    else
        averages_w=sz_w(dims_ref.subSpecs);
        rawAverages_w=1;
    end
else
    if dims_ref.averages~=0
        averages_w=sz_w(dims_ref.averages);
        rawAverages_w=averages_w;
    else
        averages_w=1;
        rawAverages_w=1;
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
if dims_ref.subSpecs~=0
    subspecs_w=sz_w(dims_ref.subSpecs);
    rawSubspecs_w=subspecs_w;
else
    subspecs_w=1;
    rawSubspecs_w=subspecs_w;
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)));
ppm=f/(Bo*42.577);

% GE data assumes the center frequency to be 4.68 ppm:
centerFreq = 4.68;
ppm=ppm + centerFreq;

t=0:dwelltime:(sz(1)-1)*dwelltime;


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
out.seq='';
out.te=GEhdr.TE;
out.tr=GEhdr.TR;
out.pointsToLeftshift=0;
out.centerFreq = centerFreq;
out.geometry = GEhdr.geometry;
temp = out.geometry.size;
out.geometry.size = [];
out.geometry.size.dim1 = temp(1);
out.geometry.size.dim2 = temp(2);
out.geometry.size.dim3 = temp(3);
if isnumeric(GEhdr.version)
    out.software = ['Rev_number ' num2str(GEhdr.version)];
else
    out.software = ['Rev_number ' GEhdr.version];
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
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end


%FOR WATER UNSUPPRESSED DATA
%FILLING IN DATA STRUCTURE
out_ref.fids=fids_ref;
out_ref.specs=specs_w;
out_ref.sz=sz_w;
out_ref.ppm=ppm;  
out_ref.t=t;    
out_ref.spectralwidth=spectralwidth;
out_ref.dwelltime=dwelltime;
out_ref.txfrq=txfrq;
out_ref.date=date;
out_ref.dims=dims_ref;
out_ref.Bo=Bo;
out_ref.averages=averages_w;
out_ref.rawAverages=rawAverages_w;
out_ref.subspecs=subspecs_w;
out_ref.rawSubspecs=rawSubspecs_w;
out_ref.seq='';
out_ref.te=GEhdr.TE;
out_ref.tr=GEhdr.TR;
out_ref.pointsToLeftshift=0;
out_ref.centerFreq = centerFreq;
out_ref.geometry = GEhdr.geometry;
temp = out_ref.geometry.size;
out_ref.geometry.size = [];
out_ref.geometry.size.dim1 = temp(1);
out_ref.geometry.size.dim2 = temp(2);
out_ref.geometry.size.dim3 = temp(3);
out_ref.software = out.software;
% Add info for niiwrite
out_ref.PatientPosition = '';
out_ref.Manufacturer = 'GE';
[~,filenamest,ext] = fileparts(filename);
out_ref.OriginalFile = [filenamest ext];
%FILLING IN THE FLAGS
out_ref.flags.writtentostruct=1;
out_ref.flags.gotparams=1;
out_ref.flags.leftshifted=0;
out_ref.flags.filtered=0;
out_ref.flags.zeropadded=0;
out_ref.flags.freqcorrected=0;
out_ref.flags.phasecorrected=0;
out_ref.flags.averaged=0;
out_ref.flags.addedrcvrs=0;
out_ref.flags.subtracted=0;
out_ref.flags.writtentotext=0;
out_ref.flags.downsampled=0;
if out_ref.dims.subSpecs==0
    out_ref.flags.isFourSteps=0;
else
    out_ref.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end

% Add info for niiwrite
out.PatientPosition = '';
out.Manufacturer = 'GE';
[~,filename,ext] = fileparts(filename);
out.OriginalFile = [filename ext];

% Sequence flags
out.flags.isUnEdited = 0;
out.flags.isMEGA = 0;
out.flags.isHERMES = 0;
out.flags.isHERCULES = 0;
out.flags.isPRIAM = 0;
out.flags.isMRSI = 0;

% Separate ecc and quant reference for CMRR sequence
if twoRefs                  
    out_w = out_ref;

    out_w.fids = out_w.fids(:,:,[3,4,7,8]);
    out_w.specs = out_w.specs(:,:,[3,4,7,8]);
    out_w.sz = size(out_w.fids);
    out_w.averages = out_w.averages/2;
    out_w.rawAverages = out_w.rawAverages/2;

    out_ref.fids = out_ref.fids(:,:,[1,2,5,6]);
    out_ref.specs = out_ref.specs(:,:,[1,2,5,6]);
    out_ref.sz = size(out_ref.fids);
    out_ref.averages = out_ref.averages/2;
    out_ref.rawAverages = out_ref.rawAverages/2;
else
    out_w = [];
end


%DONE
