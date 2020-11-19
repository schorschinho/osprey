%%read cpx
function data = read_cpx(file, border, flip_img, kspace, read_params, compression_parameter)

%--------------------------------------------------------------------------
% data = read_cpx(file, border, flip_img)
%
% read_cpx: Function for reading the whole cpx-file and writes the data in
%           a predefined 9 dimensional array
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               border      0 returns the original sized image
%                           1 makes a square image
%               flip_img    0 no flip of the image
%                           1 flips the image (like the scanner does)
%               read_params optional input. Specifies the images to be
%                           read. read_params is a struct created
%                           by the function create_read_param_struct
%
% Output:       data        A 9 dimensional array which consists all the
%                           data in the cpx file. The array has the
%                           following structure:
%                           Column: 1:resolution y, 2:resolution x 3:Stack,
%                                   4:Slice, 5:Coil, 6:Heart phase, 7:Echo,
%                                   8:Dynamics, 9:Segments, 10:Segments2
%
% Notice! For scans with body coil and surface coil in different stacks,
% e.g. SENSE reference scan, the body coil is numbered as coil no. 1 in
% stack no. 2, the surface coils are numbered as coils 2,3,4... in stack no. 1!
%------------------------------------------------------------------------


switch nargin
    case 5
        compression_parameter = [];
    case 4
        read_params = create_read_param_struct(file);
        compression_parameter = [];
    case 3
        read_params = create_read_param_struct(file);
        kspace = 0;
        compression_parameter = [];
    case 2
        flip_img = 0;
        kspace = 0;
        read_params = create_read_param_struct(file);
        compression_parameter = [];
    case 1
        flip_img = 0;
        border = 0;
        kspace = 0;
        read_params = create_read_param_struct(file);
        compression_parameter = [];
end

%Reads the header of the file
header = read_cpx_header(file,'no');
[rows,columns] = size(header);

% Calculates the number of slices, coils, etc...
stacks   = length(read_params.loca);
slices   = length(read_params.slice);
coils    = length(read_params.coil);
hps      = length(read_params.phase);
echos    = length(read_params.echo);
dynamics = length(read_params.dyn);
segments = length(read_params.seg);
segments2= length(read_params.seg2);


res_x       = header(1,9);
res_y       = header(1,10);
compression = header(1,11);
flip        = header(1,12);

offset_table_cpx=create_offset_table(header);

% defines the array for the output
if border
    res = max([res_x, res_y]);
    res_x = res;
    res_y = res;
end

if ~isempty(compression_parameter)
    data = zeros(res_x, res_y, stacks, slices, compression_parameter{1}, hps, echos, dynamics, segments, segments2,'single');
    data3 = zeros(res_x, res_y,coils);
else
    data = zeros(res_x, res_y, stacks, slices, coils, hps, echos, dynamics, segments, segments2,'single');
    data3 = zeros(res_x, res_y,coils);
end

% Define a waitbar
h = waitbar(0, 'Loading file...');
set(h,'Units','pixels')
scnsize = get(0,'ScreenSize');
set(h,'Position',[floor(scnsize(3)/2)-160,floor(scnsize(4)/2)-30,360,75]);

fid = fopen(file);

% Runs through all images in the file, reads them and writes it in the
% correct position in the output array "data"
i = 1;
total_loops = 1;
for loop = 1:2
    for st = 1:stacks
        for sl = 1:slices
            for se2 = 1:segments2
                for ph = 1:hps
                    for ec = 1:echos
                        for dy = 1:dynamics
                            for se = 1:segments
                                for co = 1:coils
                                    offset = offset_table_cpx(read_params.loca(st),read_params.slice(sl),read_params.coil(co),read_params.phase(ph),read_params.echo(ec),read_params.dyn(dy),read_params.seg(se),read_params.seg2(se2));
                                    if offset >=0
                                        if loop == 2
                                            image = read_cpx_image(file, offset, border, flip_img);
                                            if kspace
                                                image = fftshift(fft2(fftshift(image)));
                                            end
                                            data3(:,:,co) = image;
                                            waitbar(i/total_loops,h)
                                            i = i+1;
                                        else
                                            total_loops = total_loops +1;
                                        end
                                    end
                                end
                                if ~isempty(compression_parameter)
                                    %                                     data_temp= squeeze(combine_data(reshape(data3,size(data3,1),size(data3,2), 1,size(data3,3)),compression_parameter{2}));
                                    data(:,:,st,sl,:,ph,ec,dy,se, se2) = reshape(combine_data_gui(reshape(data3,size(data3,1),size(data3,2), 1,size(data3,3)),compression_parameter{2}),size(data3,1),size(data3,2),1,1,compression_parameter{1});
                                else
                                    data(:,:,st,sl,:,ph,ec,dy,se, se2) = reshape(data3,size(data3,1),size(data3,2),1,1,size(data3,3));
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

close(h);
fclose all;

end

%%create read param struct
function [v,raw_params] = create_read_param_struct(file)

dotind = findstr(file,'.');
ending = lower(file(dotind(end)+1:end));

switch ending
    case 'rec'
        parfile = [file(1:dotind),'par'];
        par = parread(parfile);
        v.slice = unique(par.ImageInformation.SliceNumber);
        v.echo = unique(par.ImageInformation.EchoNumber);
        v.dyn = unique(par.ImageInformation.DynamicScanNumber);
        v.phase =  unique(par.ImageInformation.CardiacPhaseNumber);
        v.typ = unique(par.ImageInformation.ImageTypeMr);
        raw_params = par;
    case 'par'
        par = parread(file);
        v.slice = unique(par.ImageInformation.SliceNumber);
        v.echo = unique(par.ImageInformation.EchoNumber);
        v.dyn = unique(par.ImageInformation.DynamicScanNumber);
        v.phase =  unique(par.ImageInformation.CardiacPhaseNumber);
        v.typ = unique(par.ImageInformation.ImageTypeMr);
        raw_params = par;
    case 'cpx'
        header = read_cpx_header(file,'no');
        v.loca = unique(header(:,1))+1;
        v.slice = unique(header(:,2))+1;
        v.coil = unique(header(:,3))+1;
        v.phase = unique(header(:,4))+1;
        v.echo = unique(header(:,5))+1;
        v.dyn = unique(header(:,6))+1;
        v.seg = unique(header(:,7))+1;
        v.seg2 = unique(header(:,18))+1;
        raw_params = header;
    case 'data'
        t = 'TEHROA';
        listfile = [file(1:dotind),'list'];
        list = listread(listfile);
        typ = unique(list.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = findstr(typ(i),t);
        end
        v.typ = sort(numtyp);
        v.mix = unique(list.Index.mix)+1;
        v.dyn = unique(list.Index.dyn)+1;
        v.phase = unique(list.Index.card)+1;
        v.echo = unique(list.Index.echo)+1;
        v.loca = unique(list.Index.loca)+1;
        v.coil = unique(list.Index.chan)+1;
        v.seg = unique(list.Index.extr1)+1;
        v.seg2 = unique(list.Index.extr2)+1;
        v.ky = unique(list.Index.ky)+1;
        v.slice = unique(list.Index.kz)+1;
        v.aver = unique(list.Index.aver)+1;
        raw_params = list;
    case 'list'
        t = 'TEHROA';
        list = listread(file);
        typ = unique(list.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = findstr(typ(i),t);
        end
        v.typ = sort(numtyp);
        v.mix = unique(list.Index.mix)+1;
        v.dyn = unique(list.Index.dyn)+1;
        v.phase = unique(list.Index.card)+1;
        v.echo = unique(list.Index.echo)+1;
        v.loca = unique(list.Index.loca)+1;
        v.coil = unique(list.Index.chan)+1;
        v.seg = unique(list.Index.extr1)+1;
        v.seg2 = unique(list.Index.extr2)+1;
        v.ky = unique(list.Index.ky)+1;
        v.slice = unique(list.Index.kz)+1;
        v.aver = unique(list.Index.aver)+1;
        raw_params = list;
    otherwise
        v = -1;
end
end

%%read_cpx_header
function header = read_cpx_header(file, output)

%--------------------------------------------------------------------------
% header = read_cpx_header(file)
%
% read_cpx_header: Function for reading the header of a cpx file
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               output      can be 'yes' or 'no'. Default is 'yes'
%                           Specifies if the information about the cpx file is written
%                           to the command line
%
% Output:       header      gives out the header of the cpx file as an
%                           array. The structure of this array is the
%                           following:
%                           Column: 1:Stack, 2:Slice, 3:Coil, 4:Heart phase,
%                                   5:Echo, 6:Dynamics, 7:Segments,
%                                   8:data offset, 9:Scaling factor1,
%                                   10:Scaling factor2, 11:Compression,
%                                   12:Flip, 13: Scaling factor1, 14:
%                                   scaling factor2, 15: Mix, 16: Prep Dir,
%                                   17: Sequence Number 18: Segment2,
%                                   19: Syncho Nr.
%
%------------------------------------------------------------------------
header = [];

if nargin == 1
    output = 'yes';
end

% Calculate the size of the cpx-file
fid = fopen(file);
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');

% Read in the first header for
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 111,'long');

res_x = h1(11);
res_y = h1(12);
compression = h1(14);
if ~h2(26)
    offset = h1(10);
else
    offset = h2(26);
end

matrix_data_blocks = h1(13);

% Calculates the number of images in the cpx-file
image_exist = 1; i=0;
while image_exist
    %     header_offset = (res_x * res_y * 8 /compression + offset)*i;
    header_offset = (matrix_data_blocks * 512 + offset)*i;
    fseek(fid, header_offset, 'bof');
    h1 = fread(fid, 15, 'long');
    image_exist = h1(9);
    i = i+1;
end
images = i-1;

% Defines the header:
% header Columns : 1:Stack, 2:Slice, 3:Coil, 4:Heart phase, 5:Echo, 6:Dynamics,
%                  7:Segments, 8:data offset, 9:Resolution x, 10:Resolution y,
%                  11: Compression, 12: Flip, 13:Scaling factor1, 14:Scaling factor2
%                  15: Mix, 16: Prep Dir, 17: Sequence Nr.
%                  18: Segment2, 19: Syncho Number
header = zeros(images, 19);

% Runs over all images in the file and writes out its header
for i = 0: images-1
    header_offset = (matrix_data_blocks * 512 + offset)*i;
    fseek(fid, header_offset, 'bof');
    h1 = fread(fid, 15, 'long');
    factor = fread(fid,2,'float');
    h2 = fread(fid, 111,'long');
    header(i+1,1) = h1(2);                  % Stack
    header(i+1,2) = h1(3);                  % Slice
    header(i+1,3) = h2(2);                  % Coil
    header(i+1,4) = h1(6);                  % Heart phase
    header(i+1,5) = h1(5);                  % Echo
    header(i+1,6) = h1(7);                  % Dynamics
    header(i+1,7) = h1(8);                  % Segments
    if ~h2(26)
        header(i+1,8) = h1(10);                 % Data offset
    else
        header(i+1,8) = h2(26);                 % Data offset
    end
    header(i+1,9) = h1(11);                 % Resolution x
    header(i+1,10) = h1(12);                % Resolution y
    header(i+1,11) = h1(14);                % Compression
    header(i+1,12) = h2(111);               % Flip
    header(i+1,13) = factor(1);             % Scaling factor 1
    header(i+1,14) = factor(2);             % Scaling factor 2
    header(i+1,15) = h1(1);                 % mix
    header(i+1,16) = h1(4);                 % Prep Dir
    header(i+1,17) = h1(15);                % Sequence Number
    header(i+1,18) = h2(1);                 % Segment2
    header(i+1,19) = h2(3);                 % Syncho number
    
    if h1(9) == 0
        'Header Problem!! Too many images calculated'
        break
    end
end

% Reads in the last header and checks the parameter "Complex Matrix
% Existence" if it hasn't the value 0, the file is corrupt
last_header_offset = (matrix_data_blocks * 512 + offset)*images;
fseek(fid, last_header_offset, 'bof');
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 10,'long');
if h1(9) ~= 0
    'Header Problem'
    return
end

% Prints the parameters on the screen
if strcmp(output,'yes')
    s1=sprintf('\nResolution in x-direction: %d \nResolution in y-direction: %d \nNumber of stacks: %d \nNumber of slices: %d \nNumber of coils: %d \nNumber of heart phases: %d \nNumber of echos: %d \nNumber of dynamics: %d \nNumber of segments: %d \nNumber of segments2: %d',header(1,9),header(1,10),max(header(:,1))+1,max(header(:,2))+1,max(header(:,3))+1, max(header(:,4))+1,max(header(:,5))+1,max(header(:,6))+1,max(header(:,7))+1,max(header(:,18))+1);
    disp(s1);
end
end

%%create offset table
function offset_table_cpx=create_offset_table(header)

%--------------------------------------------------------------------------
% offset_table_cpx=create_offset_table(header)
%
% offset_table_cpx: Creates an array with all image offsets in the cpx-file
%                   The offset of a certain image can be obtained by:
%                   offset_table_cpx(stack, slice, coil, heart_phase, echo, dynamic, segment, segment2)
%
% Input:		header              The header of the cpx-file. Can be obtained
%                                   with the function "read_cpx_header"
%
% Output:       offset_table_cpx    An array with the image offsets in the
%                                   cpx-file
%
%
%------------------------------------------------------------------------


[rows, columns] = size(header);

for i = 1: rows
    if header(i,8) == 0;
        offset_table_cpx(header(i,1)+1, header(i,2)+1, header(i,3)+1, header(i,4)+1, header(i,5)+1, header(i,6)+1, header(i,7)+1, header(i,18)+1) = -100;
    else
        offset_table_cpx(header(i,1)+1, header(i,2)+1, header(i,3)+1, header(i,4)+1, header(i,5)+1, header(i,6)+1, header(i,7)+1, header(i,18)+1) = header(i,8);
    end
end
offset_table_cpx(find(offset_table_cpx==0)) = -1;
offset_table_cpx(find(offset_table_cpx==-100)) = 0;
end

%%read cpx image
function data = read_cpx_image(file, offset, border, flip_img)

%--------------------------------------------------------------------------
% data = read_cpx_image(file, offset, border, flip_img)
%
% read_cpx_image : Function for reading one image, defined through the
%                  input parameters
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               offset      offset to the image from header (element (.,8))
%               border      0 returns the original sized image
%                           1 makes a square image
%               flip_img    0 no flip of the image
%                           1 flips the image (like the scanner does)
%
%               Notice! All numbers starts from 1. E.g.the first slice is
%               numbered with 1
%
% Output:       data        The requested image in a 2 dimensional array
%
%
% Notice! For scans with body coil and surface coil in different stacks,
% e.g. SENSE reference scan, the body coil is numbered as coil no. 1 in
% stack no. 2, the surface coils are numbered as coils 2,3,4... in stack no. 1!
%------------------------------------------------------------------------


%Reads the header of the requested image
fid = fopen(file);

fseek(fid, offset-512,'bof');
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 10,'long');

res_x = h1(11);
res_y = h1(12);
compression = h1(14);

%Reads the requested image
fseek(fid, offset,'bof');
switch (compression)
    case 1
        data = zeros(res_x*res_y*2,1,'single');
        data=fread(fid, res_x*res_y*2, 'float');
    case 2
        data = zeros(res_x*res_y*2,1,'single');
        data(:)=fread(fid, res_x*res_y*2,'short');
        data=factor(2)+factor(1).*data;
    case 4
        data = zeros(res_x*res_y*2,1,'single');
        data=fread(fid, res_x*res_y*2, 'int8');
        data=factor(2)+factor(1).*data;
end
data = complex(data(1:2:end),data(2:2:end));
data = reshape(data,res_x,res_y);

%Adds the border if requested
if border & (res_x ~= res_y)
    res = max([res_x, res_y]);
    data_temp = zeros(res, res);
    if res_x > res_y
        data_temp(:,floor((res - res_y)/2): res - ceil((res - res_y)/2+0.1)) = data;
    else
        data_temp(floor((res - res_x)/2): res - ceil((res - res_x)/2+0.1),:) = data;
    end
    data = data_temp;
    clear data_temp;
end

%Flips the image if requested
if flip_img
    s = size(data);
    data = data(end:-1:1,:);
    %     data=data';
    %     data = data(:,s(1)+1-(1:s(1)));
    %     data = data(s(2)+1-(1:s(2)),:);
end

fclose(fid);
end

%%load noise from raw file
function [data,info] = read_noise(filename,varargin)

% Start execution time clock and initialize DATA and INFO to empty arrays
tic;
data=[];
info=[];

% Initialize INFO structure
% Serves to fix the display order
info.filename = [];
info.loadopts = [];
info.dims = [];
info.labels = [];
info.labels_row_index_array = [];
info.label_fieldnames = [];
info.idx = [];
info.fseek_offsets = [];
info.nLabels = [];
info.nLoadedLabels = [];
info.nDataLabels = [];
info.nNormalDataLabels = [];
info.datasize = [];


% Parse the filename.
% It may be the LAB filename, RAW filename or just the filename prefix
% Instead of REGEXP, use REGEXPI which igores case
toks = regexpi(filename,'^(.*?)(\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
labname = sprintf('%s.lab',prefix);
rawname = sprintf('%s.raw',prefix);
info.filename = filename;
% Open LAB file and read all hexadecimal labels
labfid = fopen(labname,'r');
if labfid==-1,
    error( sprintf('Cannot open %s for reading', labname) );
end

% Read all hexadecimal labels
[unparsed_labels, readsize] = fread (labfid,[16 Inf], 'uint32=>uint32');
info.nLabels = size(unparsed_labels,2);
fclose(labfid);

% Parse hexadecimal labels
% Inspired by Holger Eggers' readRaw.m.  Thanks Holger! 
% See arsrcglo1.h for more details.
info.labels.DataSize.vals         = unparsed_labels(1,:);
info.labels.LeadingDummies.vals   = bitshift (bitand(unparsed_labels(2,:), (2^16-1)),  -0);
info.labels.TrailingDummies.vals  = bitshift (bitand(unparsed_labels(2,:), (2^32-1)), -16);
info.labels.SrcCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^16-1)),  -0);
info.labels.DstCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^32-1)), -16);
info.labels.SeqNum.vals           = bitshift (bitand(unparsed_labels(4,:), (2^16-1)),  -0);
info.labels.LabelType.vals        = bitshift (bitand(unparsed_labels(4,:), (2^32-1)), -16);
info.labels.ControlType.vals      = bitshift( bitand(unparsed_labels(5,:),  (2^8-1)),  -0);
info.labels.MonitoringFlag.vals   = bitshift( bitand(unparsed_labels(5,:), (2^16-1)),  -8);
info.labels.MeasurementPhase.vals = bitshift( bitand(unparsed_labels(5,:), (2^24-1)), -16);
info.labels.MeasurementSign.vals  = bitshift( bitand(unparsed_labels(5,:), (2^32-1)), -24);
info.labels.GainSetting.vals      = bitshift( bitand(unparsed_labels(6,:),  (2^8-1)),  -0);
info.labels.Spare1.vals           = bitshift( bitand(unparsed_labels(6,:), (2^16-1)),  -8);
info.labels.Spare2.vals           = bitshift (bitand(unparsed_labels(6,:), (2^32-1)), -16);
info.labels.ProgressCnt.vals      = bitshift (bitand(unparsed_labels(7,:), (2^16-1)),  -0);
info.labels.Mix.vals              = bitshift (bitand(unparsed_labels(7,:), (2^32-1)), -16);
info.labels.Dynamic.vals          = bitshift (bitand(unparsed_labels(8,:), (2^16-1)),  -0);
info.labels.CardiacPhase.vals     = bitshift (bitand(unparsed_labels(8,:), (2^32-1)), -16);
info.labels.Echo.vals             = bitshift (bitand(unparsed_labels(9,:), (2^16-1)),  -0);
info.labels.Location.vals         = bitshift (bitand(unparsed_labels(9,:), (2^32-1)), -16);
info.labels.Row.vals              = bitshift (bitand(unparsed_labels(10,:), (2^16-1)),  -0);
info.labels.ExtraAtrr.vals        = bitshift (bitand(unparsed_labels(10,:), (2^32-1)), -16);
info.labels.Measurement.vals      = bitshift (bitand(unparsed_labels(11,:), (2^16-1)),  -0);
info.labels.E1.vals               = bitshift (bitand(unparsed_labels(11,:), (2^32-1)), -16);
info.labels.E2.vals               = bitshift (bitand(unparsed_labels(12,:), (2^16-1)),  -0);
info.labels.E3.vals               = bitshift (bitand(unparsed_labels(12,:), (2^32-1)), -16);
info.labels.RfEcho.vals           = bitshift (bitand(unparsed_labels(13,:), (2^16-1)),  -0);
info.labels.GradEcho.vals         = bitshift (bitand(unparsed_labels(13,:), (2^32-1)), -16);
info.labels.EncTime.vals          = bitshift (bitand(unparsed_labels(14,:), (2^16-1)),  -0);
info.labels.RandomPhase.vals      = bitshift (bitand(unparsed_labels(14,:), (2^32-1)), -16);
info.labels.RRInterval.vals       = bitshift (bitand(unparsed_labels(15,:), (2^16-1)),  -0);
info.labels.RTopOffset.vals       = bitshift (bitand(unparsed_labels(15,:), (2^32-1)), -16);
info.labels.ChannelsActive.vals   = unparsed_labels(16,:);

clear unparsed_labels;

% Find unique values of each label field
info.label_fieldnames = fieldnames(info.labels);
for k=1:length(info.label_fieldnames),
    info.labels.(info.label_fieldnames{k}).uniq = unique( info.labels.(info.label_fieldnames{k}).vals ); 
end

% Calculate fseek offsets
info.fseek_offsets = zeros(info.nLabels,1);
info.fseek_offsets(1)=512; % add mysterious 512 byte offset to begin reading file
for k=2:info.nLabels,
    info.fseek_offsets(k) = info.fseek_offsets(k-1)+ info.labels.DataSize.vals(k-1) - info.labels.TrailingDummies.vals(k-1)  - info.labels.LeadingDummies.vals(k-1);
end
info.idx.no_data = find(info.labels.DataSize.vals==0);
info.fseek_offsets(info.idx.no_data) = -1;

% Find indices of different label control types
% See arsrcglo1.h for more details.
standard_labels = info.labels.LabelType.vals==32513;
info.idx.NORMAL_DATA         = find(info.labels.ControlType.vals== 0 & standard_labels);
info.idx.DC_OFFSET_DATA      = find(info.labels.ControlType.vals== 1 & standard_labels);
info.idx.JUNK_DATA           = find(info.labels.ControlType.vals== 2 & standard_labels);
info.idx.ECHO_PHASE_DATA     = find(info.labels.ControlType.vals== 3 & standard_labels);
info.idx.NO_DATA             = find(info.labels.ControlType.vals== 4 & standard_labels);
info.idx.NEXT_PHASE          = find(info.labels.ControlType.vals== 5 & standard_labels);
info.idx.SUSPEND             = find(info.labels.ControlType.vals== 6 & standard_labels);
info.idx.RESUME              = find(info.labels.ControlType.vals== 7 & standard_labels);
info.idx.TOTAL_END           = find(info.labels.ControlType.vals== 8 & standard_labels);
info.idx.INVALIDATION        = find(info.labels.ControlType.vals== 9 & standard_labels);
info.idx.TYPE_NR_END         = find(info.labels.ControlType.vals==10 & standard_labels);
info.idx.VALIDATION          = find(info.labels.ControlType.vals==11 & standard_labels);
info.idx.NO_OPERATION        = find(info.labels.ControlType.vals==12 & standard_labels);
info.idx.DYN_SCAN_INFO       = find(info.labels.ControlType.vals==13 & standard_labels);
info.idx.SELECTIVE_END       = find(info.labels.ControlType.vals==14 & standard_labels);
info.idx.FRC_CH_DATA         = find(info.labels.ControlType.vals==15 & standard_labels);
info.idx.FRC_NOISE_DATA      = find(info.labels.ControlType.vals==16 & standard_labels);
info.idx.REFERENCE_DATA      = find(info.labels.ControlType.vals==17 & standard_labels);
info.idx.DC_FIXED_DATA       = find(info.labels.ControlType.vals==18 & standard_labels);
info.idx.DNAVIGATOR_DATA     = find(info.labels.ControlType.vals==19 & standard_labels);
info.idx.FLUSH               = find(info.labels.ControlType.vals==20 & standard_labels);
info.idx.RECON_END           = find(info.labels.ControlType.vals==21 & standard_labels);
info.idx.IMAGE_STATUS        = find(info.labels.ControlType.vals==22 & standard_labels);
info.idx.TRACKING            = find(info.labels.ControlType.vals==23 & standard_labels);
info.idx.FLUOROSCOPY_TOGGLE  = find(info.labels.ControlType.vals==24 & standard_labels);
info.idx.REJECTED_DATA       = find(info.labels.ControlType.vals==25 & standard_labels);
info.idx.UNKNOWN27           = find(info.labels.ControlType.vals==27 & standard_labels);
info.idx.UNKNOWN28           = find(info.labels.ControlType.vals==28 & standard_labels);

% Calculate number of standard, normal data labels
info.nNormalDataLabels = length(info.idx.NORMAL_DATA);

% Dimension names
dimnames = {'coil','kx','ky','kz','E3','loc','ec','dyn','ph','row','mix','avg'};
dimfields = {'N/A','N/A','E1','E2','E3','Location','Echo','Dynamic','CardiacPhase','Row','Mix','Measurement'};

% Initialize dimension data to zero
info.dims.nCoils         = 0;
info.dims.nKx            = 0;
info.dims.nKy            = 0;
info.dims.nKz            = 0;
info.dims.nE3            = 0;
info.dims.nLocations     = 0;
info.dims.nEchoes        = 0;
info.dims.nDynamics      = 0;
info.dims.nCardiacPhases = 0;
info.dims.nRows          = 0;
info.dims.nMixes         = 0;
info.dims.nMeasurements  = 0;

% Calculate max number of active coils
maxChannelsActiveMask = 0;
for k=1:length(info.labels.ChannelsActive.uniq),
    maxChannelsActiveMask = bitor(maxChannelsActiveMask,info.labels.ChannelsActive.uniq(k));
end
while maxChannelsActiveMask > 0
    if bitand(maxChannelsActiveMask, 1),
        info.dims.nCoils = info.dims.nCoils + 1;
    end
    maxChannelsActiveMask = bitshift (maxChannelsActiveMask, -1);
end

% Calculate dimensions of normal data
info.dims.nKx            = max(info.labels.DataSize.vals(info.idx.NORMAL_DATA)) / info.dims.nCoils / 2 / 2;
info.dims.nKy            = length(unique(info.labels.E1.vals(info.idx.NORMAL_DATA)));
info.dims.nKz            = length(unique(info.labels.E2.vals(info.idx.NORMAL_DATA)));
info.dims.nE3            = length(unique(info.labels.E3.vals(info.idx.NORMAL_DATA)));
info.dims.nLocations     = length(unique(info.labels.Location.vals(info.idx.NORMAL_DATA)));
info.dims.nEchoes        = length(unique(info.labels.Echo.vals(info.idx.NORMAL_DATA)));
info.dims.nDynamics      = length(unique(info.labels.Dynamic.vals(info.idx.NORMAL_DATA)));
info.dims.nCardiacPhases = length(unique(info.labels.CardiacPhase.vals(info.idx.NORMAL_DATA)));
info.dims.nRows          = length(unique(info.labels.Row.vals(info.idx.NORMAL_DATA)));
info.dims.nMixes         = length(unique(info.labels.Mix.vals(info.idx.NORMAL_DATA)));
info.dims.nMeasurements  = length(unique(info.labels.Measurement.vals(info.idx.NORMAL_DATA)));

% With known possible dimension names, the load options can now be parsed
p = inputParser;
p.StructExpand = true;
p.CaseSensitive = true;
p.KeepUnmatched = false; % throw an error for unmatched inputs
p.addRequired('filename', @ischar);
for k=1:length(dimnames),
    p.addParamValue(dimnames{k}, [], @isnumeric);
end
p.addParamValue('verbose', false, @islogical);
p.addParamValue('savememory', true, @islogical);
p.parse(filename, varargin{:});

% Return loadopts structure inside INFO structure
% remove filename field - it is passed as the first required argument
info.loadopts = rmfield(p.Results,'filename');

% Find the unique set of values for each dimension name
info.dims.coil = [1:info.dims.nCoils];
info.dims.kx   = [1:info.dims.nKx];
for k=3:length(dimnames), % skip coil and kx
    info.dims.(dimnames{k}) = unique(info.labels.(dimfields{k}).vals(info.idx.NORMAL_DATA));
end

% Find intersection of available dimensions with LOADOPTS dimensions
for k=1:length(dimnames),
    if ~isempty(info.loadopts.(dimnames{k})),
        info.dims.(dimnames{k}) = intersect_a_with_b(info.loadopts.(dimnames{k}),info.dims.(dimnames{k}));
    end
end

% Calculate data size
datasize = []; 
for k=1:length(dimnames),
    datasize = [datasize length(info.dims.(dimnames{k}))];
end
info.datasize = datasize;

% throw error if any dimension size is zero
if any(info.datasize==0),
    zero_length_str = sprintf(' ''%s'' ', dimnames{find(info.datasize==0)});
    error('size of selected data to load has zero length along dimension(s): %s', zero_length_str);
end

% Skip data loading if only one output argument is provided, return INFO
if nargout==1,
    info.labels_row_index_array = [1:size(info.labels,1)];
    data=info;
    return;
end

% Create array to hold label row numbers for loaded data
% skip the coil and kx dimensions
info.labels_row_index_array = zeros(datasize(3:end));

% Pre-allocate DATA array
if info.loadopts.savememory==true,
    data = zeros(info.datasize,'single');
else
    data = zeros(info.datasize);
end

% Read RAW data for selected dimension ranges
fidraw = fopen(rawname,'r','ieee-le');
if fidraw<0,
    error(sprintf('cannot open RAW file: %s', rawname));
end
info.nLoadedLabels=0;

raw_data_fread_size = double(info.dims.nCoils * info.dims.nKx * 2);
rawdata_2d = complex(zeros(info.dims.nCoils,info.dims.nKx),zeros(info.dims.nCoils,info.dims.nKx));

% Read FRC noise data
if(~isempty(info.idx.FRC_NOISE_DATA))
    
    for n=1:length(info.idx.FRC_NOISE_DATA)
        
        frc_noise_idx = info.idx.FRC_NOISE_DATA(n);
        
        num = double(info.labels.ChannelsActive.vals(frc_noise_idx));
        [f,e]=log2(num);
        ncoils = sum(rem(floor(num*pow2(1-e:0)),2));
        
        frc_noise_samples_per_coil = info.labels.DataSize.vals(frc_noise_idx) / 2 / 2 / ncoils;
    byte_offset = info.fseek_offsets(frc_noise_idx);
    status = fseek(fidraw, byte_offset, 'bof');
    rawdata_1d = double(fread(fidraw, double(info.labels.DataSize.vals(frc_noise_idx)/2) , 'int16'));
    info.FRC_NOISE_DATA(1:ncoils,:,n) = permute(reshape(rawdata_1d(1:2:end) + 1i*rawdata_1d(2:2:end), frc_noise_samples_per_coil, ncoils),[2 1]);
    end
end
fclose(fidraw);

% Calculate total raw data blobs
size_data = size(data);
max_img_dims = size_data(3:end);
info.nDataLabels = prod(max_img_dims);

% If VERBOSE, display execution information
if info.loadopts.verbose==true,
    disp( sprintf('Loaded %d of %d available normal data labels', info.nLoadedLabels, info.nNormalDataLabels) );
    tmpstr = '';
    for k=1:length(dimnames),
        tmpstr = sprintf('%s, # %s: %d', tmpstr, dimnames{k}, length(info.dims.(dimnames{k})) );
    end
    disp( sprintf('Data contains %d raw labels - %s', info.nDataLabels, tmpstr(3:end)) );
    disp( sprintf('Total execution time = %.3f seconds', toc) );
end

% Find intersection of vector a with vector b without sorting 
function c = intersect_a_with_b(a,b)
c = a;
% work backwards in order to use [] assignment
for k=length(a):-1:1,
    if isempty(find(a(k)==b)),
        c(k)=[]; 
    end
end

% force c to be a row vector for easier display
c = c(:).';


end
end