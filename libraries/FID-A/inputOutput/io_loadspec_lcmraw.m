function out = io_loadspec_lcmraw(rawFile)
%% LCMparam = io_loadspec_lcmraw(rawFile)
%   This function reads a LCModel-formatted .RAW file.
%
%   USAGE:
%      LCMparam     = osp_readlcm_raw(rawFile)
%
%   INPUTS:
%      rawFile      = LCModel .RAW data file
%
%   OUTPUTS:
%      % out        = Input dataset in FID-A structure format.
%
%   AUTHORS:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-08-16)
%       goeltzs1@jhmi.edu
%   

% Throw error if no input is provided
if nargin < 1
    error('No raw file provided to read.');
end

% Open the provided file, create empty cell to store the parsed entries
fid = fopen(rawFile);
C = {};

% The first line of the .RAW header is always '$SEQPAR', the last line always '$END'
firstLine = '$SEQPAR';
lastLine  = '$END';

% Get lines and continue until we get to the last line indicator
tline = fgets(fid);
while isempty(strfind(tline, lastLine))
    tline = fgets(fid);
    
    % Parse every line between the start and the end line
    if (isempty(strfind(tline ,firstLine)) + isempty(strfind(tline, lastLine)) == 2)
        C{end+1}= strtrim(strsplit(tline, {'=',','}));
    end
    
end

% The first line of the NMID namelist is always '$NMID', the last line always '$END'
firstLine = '$NMID';
lastLine  = '$END';

% Get lines and continue until we get to the last line indicator
tline = fgets(fid);
while isempty(strfind(tline, lastLine))
    tline = fgets(fid);
    
    % Parse every line between the start and the end line
    if (isempty(strfind(tline ,firstLine)) + isempty(strfind(tline, lastLine)) == 2)
        C{end+1}= strtrim(strsplit(tline, {'=',','}));
    end
    
end

% Loop over everything that has been parsed
LCMparam = struct;
for rr = 1:length(C)
    % Most entries should be pairs, but there are cases where control
    % parameters have been entered on the same line and separated by
    % commas, or parameters like 'title' where the value may contain the
    % delimiters '=' and ','
    if length(C{rr}) == 2
        LCMparam = parseControlFileLine(LCMparam, C{rr});
    else
        % If the title field is split, join back together
        if strcmpi(C{rr}{1}, 'title')
            title = C{rr}{2};
            for kk = 3:length(C{rr})
                title = [title, C{rr}{kk}];
            end
            LCMparam.title = title;
            
        else
            
            % If the number of entries in this line is even, proceed, otherwise
            % throw an error
            if mod(length(C{rr}),2) == 0
                % Evaluate pairwise
                for pp = 1:length(C{rr})/2
                    P{1} = C{rr}{2*pp-1};
                    P{2} = C{rr}{2*pp};
                    LCMparam = parseControlFileLine(LCMparam, P);
                end
            else
                error('Invalid control file line: %s', C{rr});
            end
        
        end
        
    end   

end





% Now get the FID
tline=fgets(fid);
fileEnd = false;
linenum=1;
RF=[];
% If the line is empty skip it
while tline~=-1
    %dataline=line(1:semicol_index-2);
    [A,count, errmsg, nextindex] = sscanf(tline, '%f', inf);
    % If read failed, output the error     
    if ~isempty(errmsg);
       fclose(fid);
       error('READLCMRAW_BASIS failed with read error: %s', errmsg);
    end
    % Store the read values into rf array
    RF(linenum) = A(1)+1i*A(2);
    linenum = linenum + 1;
    tline=fgets(fid);
end

    out.fids=RF';
    out.specs=fftshift(fft(out.fids));
    vectorsize = size(out.fids);
    sz=[vectorsize 1];
    dwelltime = str2double(LCMparam.DELTAT); % dwelltime [s]
    spectralwidth = 1/dwelltime; % bandwidth [Hz]
    hzpppm = str2double(LCMparam.HZPPPM);
    txfrq=hzpppm*1e6; % transmitter frequency [Hz]
    Bo = txfrq/42577000; % B0 [T]
    
    f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
    ppm=f/(Bo*42.577);
    
    centerFreq = 4.68;
    ppm=ppm + centerFreq;
    t=[dwelltime:dwelltime:vectorsize*dwelltime];
    
    averages = 1;
    rawAverages = 1;
    subspecs = 1;
    rawSubspecs = 1;
    dims.t = 1;
    dims.averages = 0;
    dims.coils = 0;
    dims.extras = 0;
    dims.subSpecs = 0;
    
    geometry.size.ap = 0; % voxel size in AP direction [mm]
    geometry.size.lr = 0; % voxel size in LR direction [mm]
    geometry.size.cc = 0; % voxel size in CC direction [mm]
    geometry.pos.ap  = 0; % voxel center offset in AP direction [mm]
    geometry.pos.lr  = 0; % voxel center offset in LR direction [mm]
    geometry.pos.cc  = 0; % voxel center offset in CC direction [mm]
    geometry.rot.ap  = 0; % angulation around AP axis [deg]
    geometry.rot.lr  = 0; % angulation around LR axis [deg]
    geometry.rot.cc  = 0; % angulation around CC axis [deg]

    %FILLING IN DATA STRUCTURE
    out.sz=sz;
    out.ppm=ppm;
    out.t=t;
    out.spectralwidth=spectralwidth;
    out.dwelltime=dwelltime;
    out.txfrq=txfrq;
    out.dims=dims;
    out.Bo=Bo;
    out.averages=averages;
    out.rawAverages=rawAverages;
    out.subspecs=subspecs;
    out.rawSubspecs=rawSubspecs;
    out.seq=LCMparam.SEQ;
    out.te=str2double(LCMparam.ECHOT);
    out.tr=0; % for lack of better input
    out.pointsToLeftshift=0;
    out.centerFreq = centerFreq;
    out.geometry = geometry;
    
    %FILLING IN THE FLAGS
    out.flags.writtentostruct=1;
    out.flags.gotparams=1;
    out.flags.leftshifted=0;
    out.flags.filtered=0;
    out.flags.zeropadded=0;
    out.flags.freqcorrected=0;
    out.flags.phasecorrected=0;
    out.flags.averaged=1;
    out.flags.addedrcvrs=1;
    out.flags.subtracted=0;
    out.flags.writtentotext=0;
    out.flags.downsampled=0;
    if out.dims.subSpecs==0
        out.flags.isISIS=0;
    else
        out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
    end

    
    
fclose(fid);

end

function LCMparam = parseControlFileLine(LCMparam, keyValuePair)
    % Test whether the input key has parentheses in it, for example chomit
    key = keyValuePair{1};
    expr = '(\w+)';
    [tokens,matches] = regexp(key,expr,'tokens','match');
    % If this regular expression search returns two tokens, it means that
    % the control parameter name is a vector. We'll save the entries as a
    % cell.
    if length(tokens) == 1
        LCMparam.(key) = keyValuePair{2};
    elseif length(tokens) == 2
        LCMparam.(tokens{1}{1}){str2num(tokens{2}{1})} = keyValuePair{2};
    end
    
end
