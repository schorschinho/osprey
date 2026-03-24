function kspace = loadRawKspace(listordatafile, varargin)
% LOADRAWKSPACEPHILIPS Loads the Philips *.list and *.data files that can
% be exported from Philips MRI systems.
%
% [KSPACE] = LOADRAWKSPACEPHILIPS(LISTORDATAFILE)
% [KSPACE] = LOADRAWKSPACEPHILIPS(LISTORDATAFILE,COIL_COMBINATIONS)
%
% KSPACE            Struct containing all kspace parameters read from the 
%                   header file (.list) and binary data file (.data).
%                   More information see below.
% LISTORDATAFILE    valid .list or .data file. Note that the data and list file should be 
%                   located in the same folder.
% COIL_COMBINATIONS Cell array with doubles (or arrays of doubles) defining
%                   the coil elements that will be merged (linear
%                   combination:
%                   sum(coil_elements_data)/numel(coilelements)
%
% kspace.typ          type kline (STD,REJ,NOI)
% kspace.mix          mixed sequence number
% kspace.dyn          dynamic scan number
% kspace.card         cardiac phase number
% kspace.echo         echo number
% kspace.loca         location number
% kspace.chan         synco channel number
% kspace.extr1        extra attribute 1
% kspace.extr2        extra attribute 2
% kspace.[ky,kz]      k-space location in 1st and 2nd preparation direction
%                      (image data)
% kspace.[kx,ky,kz]   k-space location in 1st, 2nd and 3rd preparation direction
%                      (spectroscopy data)
% kspace.n_a_         ?
% kspace.aver         sequence number of this signal average
% kspace.sign         sign of measurement gradient used for this data
%                      vector (1=pos, -1=neg)
% kspace.rf           sequence number of this rf echo
%                       (only for TSE,TFE,GraSE)
% kspace.grad         sequence number of this gradient echo
%                       (only for EPI, GraSE)
% kspace.enc          encoding time (only for EPI, GraSE)
% kspace.rtop         R-top offset in ms
% kspace.rr           RR interval in ms
% kspace.size         data vector size in bytes (1 complex elements = 2
%                       single floats = 8 bytes)
% kspace.offset       data vector offset in bytes (first data vector starts
%                       at offset 0)
% kspace.complexdata  calculated complex kspace data (kx complex numbers at location ky,kz)
%
% kspace.kspace_properties  Additional values found in the header file, such
%                            as kspace offsets.
% kspace.kspace_properties.coil_combinations: cell array defining the
%    sorted input coil_combinations cell array (.coil_combinations{1,:}) 
%    and the final used cell array with coil ids (.coil_combinations{2,:})
%
% For more information on the header, please refer to the *.list file using an arbitrary text editor.
%
%
% Example:
% 
% kspace = loadRawKspacePhilips('C:/raw_004.list')
% kspace = loadRawKspacePhilips('C:/raw_004.list',{[0:3],[4:5],7,8})
%            with coil_combinations {[0:3],[4:5],7,8} a cell array defining 
%            the coils to combine.
%
% --- Author information
% Wouter Potters
% Academic Medical Center, Amsterdam, The Netherlands
% w.v.potters@amc.nl
% Date: 12-August-2014


if ischar(listordatafile)
    switch listordatafile(end-4:end)
        case '.list'
            listfile = listordatafile;
            datafile = [listordatafile(1:end-5) '.data'];
        case '.data'
            datafile = listordatafile;
            listfile = [listordatafile(1:end-5) '.list'];
        otherwise
            error('MATLAB:loadRawKspace:wrong_inputfile','Please provide a valid list or data file');
    end
end

coil_combinations = []; % no coil combinations = save all coils seperately.
                        % default behaviour
                            
if nargin > 1 % coil combinations requested
    if ischar(varargin{1})
        
    else 
    coil_combinations = varargin{1}; % cell input with coil combinations
                                     % e.g. {[0,1,2,3],[4,5,6],[7,8,9,10]}
                                     % will combine coils [0,1,2,3] and
                                     % [4,5,6] and [7,8,9,10] into one
                                     % channel.
                                     % note that the combined result will
                                     % be saved in the first coil of the
                                     % array (0, 4 and 7 in this case)
     if ~iscell(coil_combinations) || isempty(coil_combinations)
         error('Non-empty cell type expected for coil_combinations');
     elseif ~all(cellfun(@(c) isnumeric(c) & ~isempty(c) & (all(round(c) == c)),coil_combinations))
         error('Cells in coil_combinations can only contain numeric integer values')
     end
    end

end
    

% load the listfile
if exist(listfile,'file')==2
    kspace = loadListFile(listfile);
else
    error(['List file (' listfile ') does not exist'])
end
if exist(datafile,'file')==2
    kspace = loadDataFile(datafile,kspace,coil_combinations); % coil_combinations is empty if not defined
else
    warning(['Data file (' listfile ') not found. Loading the raw data was skipped.'])
end




function kspace = loadDataFile(datafile,kspace,coil_combinations)
% kspace = loadDataFile(datafile,kspace)
kspace.complexdata = []; %create extra variable.
nroflocations = numel(kspace.typ);

% read NOI_binary_data
try 
    fid = fopen(datafile);
catch %#ok<*CTCH>
    error('MATLAB:loadDataFile:invalidFile','Data file is invalid.');
end
if ( fid == -1 )
    error('MATLAB:loadDataFile:invalidFile','Data file is invalid.');
end
if isempty(coil_combinations)
     % default usage, just load all coils, no decrease in data size
    for cur = 1:nroflocations
        kspace.complexdata{cur,1} = [1 1i] * fread(fid,[2 kspace.size(cur)/8],'single',0,'ieee-le');
    end
    
elseif iscell(coil_combinations) && ~isempty(coil_combinations)
    % This piece of code assumes that the same channels are always used for
    % acquisitions of the same data and subsequently combines some of the
    % coil elements into a single (virtual) coil element
    
    % This assumption is then used to quickly read the data of all channels
    % and immediately the to-be-combined channels are combined while reading 
    % the data. This saves precious memory (#channels/#groups)
    if ~all(cellfun(@(c) size(coil_combinations,1) == 1,coil_combinations))
        error('The coil combination input should consist of matrices of size 1xN')
    end
    
    % the coil combination cell elements are sorted ascending.
    coil_combinations = cellfun(@(c) sort(c,'ascend'),coil_combinations,'uniformoutput',false);

    
    coil_combis = cat(2,coil_combinations{:});
    if ( length(unique(coil_combis)) ~= length(coil_combis) )
        error('The coil combinations should not contain duplicate coils across combinations')
    end
    
    if ( length(unique(coil_combis)) ~= length(unique(kspace.chan)) )
        error(['Not all coils were used! Missing coils: ' num2str(kspace.chan(~ismember(unique(kspace.chan),coil_combis))') '.'])
    end
    
    cur = 1; % start at the beginning
    prev_mat_size=0;
    nr_of_channels = length(unique(kspace.chan));
	index_of_delete_these_kspace_properties = find(strcmp(kspace.typ,'REJ'))'; % delete rejected lines later... still need them while loading data!
    warning(sprintf(['All rejected lines (typ=REJ, ' num2str(numel(index_of_delete_these_kspace_properties)) ' of ' num2str(numel(kspace.typ)) ' lines) were ignored while loading and combining the data. \nReason: REJ lines only contain 1 coil element.'])); %#ok<SPWRN>
    t1 = tic;
    list_of_rejected_klines = find(strcmp(kspace.typ,'REJ'))';
    while cur < nroflocations % continue while points are available
        if any(list_of_rejected_klines == cur)
            fread(fid,[2 kspace.size(cur)/8],'single',0,'ieee-le');
            % note that REJECTED (REJ type)
            % lines only contain 1 coil element
            % and are therefore ignored.
            cur = cur+1; 
            continue; % skip this line...
        end
        cur_pos_one = cur;
        mat_size = kspace.size(cur)/8;

        % only create new all_coils-matrix if lenght of freq.readout changed.
        if mat_size ~= prev_mat_size 
            kspace_all_coils = zeros(nr_of_channels,kspace.size(cur)/8); % preallocate matrix for current kspace line
        end

        % now read all the coil elements.
        for ichan = 1:nr_of_channels
            kspace_all_coils(ichan,:) = [1 1i] * fread(fid,[2 kspace.size(cur)/8],'single',0,'ieee-le');
            cur = cur + 1; % next line
        end
        
        % and average over the defined coil_combinations
        for icomb = 1:length(coil_combinations)
            kspace.complexdata{cur_pos_one+coil_combinations{icomb}(1),1} = sum(kspace_all_coils( coil_combinations{icomb}+1, : ), 1)./numel(coil_combinations{icomb}); % store linear combination of coil elements 
                                                                                                                                                                       % at position of first coil element
            index_of_delete_these_kspace_properties = [index_of_delete_these_kspace_properties cur_pos_one + coil_combinations{icomb}(2:end)]; %#ok<AGROW>
        end
        prev_mat_size = mat_size;
        if rem(cur,250*length(coil_combinations)) == 0
            time_now = toc(t1);
            disp(['#' num2str(cur) ' time_remaining: ' num2str(nroflocations * (time_now / cur) - time_now) 's (' num2str((nroflocations * (time_now / cur) - time_now)/60) ' min)']);
        end
    end

    % make sure that the length of complexdata is sufficient for deletion
    % later on...
    if (length(kspace.complexdata) < nroflocations)
        kspace.complexdata{nroflocations} = [];
    end
    
    % Delete unused items using index_of_delete_these_kspace_properties
    for ifieldname = fieldnames(kspace)'
       if ~strcmp(ifieldname{1},'kspace_properties')
            kspace.(ifieldname{1})(index_of_delete_these_kspace_properties) = [];
       end
    end
    
    % save used coil_combinations
    kspace.kspace_properties.coil_combinations = coil_combinations;
    
    % put coil ids at 0 to #coil combinations in
    % coil_combinations(2,:) to log coil id changes
    new_ids = 0:length(coil_combinations)-1;
    old_ids = sort(cellfun(@(c) c(1),coil_combinations));
    [~, new_index_newids] = sort(old_ids);
    new_ids = new_ids(new_index_newids);
    kspace.kspace_properties.coil_combinations(2,:) = num2cell(new_ids);
    
    % change the coil ids in the kspace struct
    dummy = kspace.chan;
    for ichan_rename = 1:length(old_ids)
        dummy(kspace.chan == old_ids(ichan_rename)) = new_ids(ichan_rename);
    end
    kspace.chan = dummy;

else
    error('Unexpected situation - empty cell object encountered for coil_combinations')
end
fgetl(fid); %go one line further to reach end-of-file

if feof(fid)
    fclose(fid);
else
    n = 0;
    fprintf('=============================================\n')
    while ~feof(fid)
        fgetl(fid)
        n=n+1;
    end
    warning('MATLAB:loadDataFile:DidNotReachEnd',['End of file finally reached after ' num2str(n) ' extra lines. Please make sure all data is available.']);
    fprintf('=============================================\n')
end

function kspace = loadListFile(listfile)
%kspace = loadListfile(filepath)
headers = 'typ mix   dyn   card  echo  loca  chan  extr1 extr2 kx    ky    kz  aver  sign  rf    grad  enc   rtop  rr    size   offset';
headers = textscan(headers,'%s');
headers = headers{1}(:)';
[nrofheaderlines, kspace_properties] = getNrofheaderlines(listfile);

try 
    fid = fopen(listfile);
catch
    error('MATLAB:loadListFile:invalidFile','List file is invalid.');
end


list = textscan(fid,'%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',inf,'HeaderLines',nrofheaderlines,'MultipleDelimsAsOne',1);
fclose(fid);

% remove
while strcmp(list{:,1}(end),'#')
    for i = length(list):-1:1
        try
            list{:,i}(length(list{:,1})) = [];
        catch
            %ignore error
        end
    end
end

kspace = struct(); %create empty struct
%fill struct with information on klines
for i = 1:length(headers)
    eval(['kspace.' strrep(headers{i},'.','_') '= list{:,' num2str(i) '};'])
end
kspace.complexdata = [];
kspace.kspace_properties = kspace_properties;

function [nrofheaderlines,kspace_properties] = getNrofheaderlines(listfile)
fid = fopen(listfile);
if (fid == -1)
    error('MATLAB:loadRawKspace:wrong_inputfile','Please provide a valid list file');
end
nrofheaderlines = 0;
currentline = '#';
kspace_properties = struct();
while any(strcmp(currentline(1),{'#','.'}))
   currentline = fgetl(fid);
   nrofheaderlines = nrofheaderlines + 1;
   if strcmp(currentline(1),'.')
       % save property;
       tmpi=strfind(currentline, ':');
       C = textscan(currentline(19:tmpi-1), '%s');
       C{1}{1}=strrep(C{1}{1},'-','_');
       tmpname = C{1}{1};
       if numel(C{1}) > 2
           for tmpj = 2:numel(C{1})
               tmpname = [tmpname, '_', C{1}{tmpj}]; % fix by susan to include spaces.
           end
       end
       fieldname = tmpname;
       if isfield(kspace_properties,fieldname) % for multipliple occasions...
           kspace_properties.(fieldname)(end+1,:) = str2num(currentline(tmpi+1:end));
       else
           kspace_properties.(fieldname) = str2num(currentline(tmpi+1:end));
       end
%        C = textscan(currentline,'. %f %f %f %s : %f %f',1);
%        fieldname = strrep(C{4}{1},'-','_');
%        kspace_properties.(fieldname) = [C{5} C{6}];
   end
end

nrofheaderlines = nrofheaderlines - 1;
fclose(fid);