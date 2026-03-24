function TWIXDeIdentify(fnames)
%% TWIXDeIdentify(fnames)
%   Reads Siemens TWIX files (*.dat) and removes participant information. New
%   de-identified TWIX files are then output, with filenames appended with
%   '_noID'. The original files are not overwritten.
%   
%   If TWIXDeIdentify has already been run, and de-identified files are
%   found within the directory, the user will be asked whether these files
%   should be overwritten.
%
%   NOTE: The user must make sure that filenames themselves do not contain
%   information that can personally identify participants. This function
%   will only de-identify the content of the TWIX file.
%
%   NOTE: This code has been tested with TWIX files from
%       VB17A
%       VB20P
%       VD13D
%       VE11A
%       VE11C
%
%   However, for your particular configuration/software version, it can not be
%   guaranteed that the code is running without error, or catching all
%   files that are to be de-identified. Please check back if you are unsure
%   whether de-identification has worked for you.
%
%   Usage:
%       TWIXDeIdentify, by itself, de-identifies all TWIX files found
%       within the current directory.
%
%       TWIXDeIdentify(fnames) de-identifies TWIX files listed in the
%       cell array fnames.
%
%   Example:
%       c = {'MRS_01.dat', 'MRS_02.dat'};
%       TWIXDeIdentify(c);
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2016-08-25)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
%       This code is a modification of previous work:
%       PhilipsDeIdentify.m 
%       (Dr. Mark Mikkelsen, Johns Hopkins University)
%
%       This code is heavily relying on TWIX parsing code bits from:
%       mapVBVD.m
%       read_twix_hdr.m
%       loop_mdh_read.m
%       (Dr. Philipp Ehses, University of Tuebingen)
%
%   History:
%       2016-08-25: First version of the code.
%       2016-08-26: Added various fields to be de-identified.
%       2017-02-02: Fixed errors for certain header formats.
%                   Added various fields to be de-identified.
%       2017-02-07: Added further fields to be de-identified.
%       2018-09-14: Added further fields to be de-identified.
%       2018-09-25: Minor bug fix.

if nargin < 1 % De-identify all DAT files in current directory

    flist = dir('*.dat');
    flist = flist(cellfun(@isempty, strfind({flist.name}, '._'))); %#ok<*STRCLFH>
    for ii = 1:length(flist)
        fnames(ii) = cellstr(flist(ii).name);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
else % De-identify DAT files user has listed in fnames

    % Check if filenames include a .dat extension
    for ii = 1:length(fnames)
        ext = fnames{ii}(end-3:end);
        assert(strcmpi(ext, '.dat'), ...
            ['The filename ' fnames{ii} ' (' num2str(ii) ')' ' in ' inputname(1) ...
            ' does not include a .dat extension.']);
    end

    % Check if files can be found
    for ii = 1:length(fnames)
        assert(any(exist(fnames{ii}, 'file')), ...
            ['The file ' fnames{ii} ' (' num2str(ii) ')' ' cannot be found.' ...
            ' Check spelling of filenames in ' inputname(1) ...
            ' (TWIX files must include a .dat extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end

    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
end

% Loop over all input DAT files, make copies, and edit these copies.
for ii = 1:length(fnames)

    % Make a copy of the DAT file that is to be worked on.
    fnames_noid = cell(size(fnames));
    fnames_noid{ii} = [fnames{ii}(1:end-4) '_noID' fnames{ii}(end-3:end)];
    copyfile(fnames{ii}, fnames_noid{ii});
    % Open this new file.
    fid = fopen(fnames_noid{ii}, 'r+','l','US-ASCII');

    % *** TAKEN FROM MAPVBVD.M ***
    % start of actual measurement data (sans header)
    fseek(fid,0,'bof');

    firstInt  = fread(fid,1,'uint32');
    secondInt = fread(fid,1,'uint32');

    % lazy software version check (VB or VD?)
    if and(firstInt < 10000, secondInt <= 64)
        version = 'vd';
        disp('Software version: VD (!?)');

        % number of different scans in file stored in 2nd in
        NScans = secondInt;
        measID = fread(fid,1,'uint32');
        fileID = fread(fid,1,'uint32');
        % measOffset: points to beginning of header, usually at 10240 bytes
        measOffset = fread(fid,1,'uint64');
        measLength = fread(fid,1,'uint64');
        fseek(fid,measOffset,'bof');
        hdrLength  = fread(fid,1,'uint32');

    else
        % in VB versions, the first 4 bytes indicate the beginning of the
        % raw data part of the file
        version  = 'vb';
        disp('Software version: VB (!?)');
        measOffset = 0;
        hdrLength  = firstInt;
        NScans     = 1; % VB does not support multiple scans in one file
    end
    
    datStart = measOffset + hdrLength;
    
    %SRY read data correction factors
    % do this for all VB datasets, so that the factors are available later
    % in the image_obj if the user chooses to set the correction flag
    if (strcmp(version, 'vb')) % not implemented/tested for vd, yet
        frewind(fid);
        while ( (ftell(fid) < datStart) && ~exist('rawfactors', 'var'))
            line = fgetl(fid);
            %find the section of the protocol
            %note: the factors are also available in <ParamArray."CoilSelects">
            %along with element name and FFT scale, but the parsing is
            %significantly more difficult
            if (~isempty(strfind(line, '<ParamArray."axRawDataCorrectionFactor">')))
                while (ftell(fid) < datStart)
                    line = fgetl(fid);
                    %find the line with correction factors
                    %the factors are on the first line that begins with this
                    %pattern
                    if (~isempty(strfind(line, '{ {  { ')))
                        line = strrep(line, '}  { } ', '0.0');
                        line = strrep(line, '{', '');
                        line = strrep(line, '}', '');
                        rawfactors = textscan(line, '%f');
                        rawfactors = rawfactors{1}; %textscan returns a 1x1 cell array
                        % this does not work in this location because the MDHs
                        % have not been parsed yet
                        %                    if (length(rawfactors) ~= 2*max(image_obj.NCha))
                        %                       error('Number of raw factors (%f) does not equal channel count (%d)', length(rawfactors)/2, image_obj.NCha);
                        %                    end;
                        if (mod(length(rawfactors),2) ~= 0)
                            error('Error reading rawfactors');
                        end
                        %note the transpose, this makes the vector
                        %multiplication during the read easier
                        arg.rawDataCorrectionFactors = rawfactors(1:2:end).' + 1i*rawfactors(2:2:end).';
                        break;
                    end
                end
            end
        end
        disp('Read raw data correction factors');
    end

    cPos = measOffset;
    
    % Loop over all scans within the TWIX file (only VD and higher!)
    for s=1:NScans
        
        fseek(fid,cPos,'bof');
        hdr_len = fread(fid, 1,'uint32');
        
        % Call routine to read in the header data and de-identify them.
        change_twix_hdr(fid);
        
        % Jump to first mdh
        cPos = cPos + hdr_len;
        fseek( fid, cPos, 'bof' );
        
        % Find out where this scan ends.
        [~, filePos, ~] = loop_mdh_read( fid, version );
        
        % Jump to end
        cPos = filePos( end );
        filePos( end ) = [];
    end
    fclose(fid);
    
    % Perform sanity check on file size (should be unchanged):
    bytes_old = cell(size(fnames));
    bytes_noid = cell(size(fnames));
    bytes_old{ii} = getfield(dir(fnames{ii}), 'bytes');
    bytes_noid{ii} = getfield(dir(fnames_noid{ii}), 'bytes');
    if bytes_old{ii} ~= bytes_noid{ii}
        msg = sprintf('*** WARNING! *** Filesizes of %s (%i) and %s (%i) do not agree! Check for errors!',fnames{ii},fnames_noid{ii},bytes_old{ii},bytes_noid{ii});
    else
        msg = sprintf('*** SUCCESSFULLY DEIDENTIFIED! *** Filesizes of %s (%i) and %s (%i) agree!',fnames{ii},bytes_old{ii},fnames_noid{ii},bytes_noid{ii});
    end
    disp(msg);
end


end

function [mdh_blob, filePos, isEOF] = loop_mdh_read( fid, version )
%% This is a modified version of loop_mdh_read by Philipp Ehses (University of Tübingen)
% This slightly shortened code just rushes through the mdh bits
% of the TWIX file to determine the bit where they end, i.e. where the next scan
% begins with a new header of its own to be de-identified.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2016-08-25)
%   
%   Credits:    
%       This code is heavily relying on TWIX parsing code bits from:
%       loop_mdh_read.m
%       (Dr. Philipp Ehses, University of Tuebingen)

    switch version
        case 'vb'
            isVD    = false;
            byteMDH = 128;
        case 'vd'
            isVD    = true;
            byteMDH = 184;
            szScanHeader    = 192; % [bytes]
            szChannelHeader =  32; % [bytes]
        otherwise
            % arbitrary assumptions:
            isVD    = false;
            byteMDH = 128;
            warning( [mfilename() ':UnknownVer'], 'Software version "%s" is not supported.', version );
    end

    cPos            = ftell(fid);
    n_acq           = 0;
    allocSize       = 4096;
    ulDMALength     = byteMDH;
    isEOF           = false;

    mdh_blob = zeros( byteMDH, 0, 'uint8' );
    szBlob   = size( mdh_blob, 2 );
    filePos  = zeros(0, 1, class(cPos));  % avoid bug in Matlab 2013b: https://scivision.co/matlab-fseek-bug-with-uint64-offset/

    % get file size
    fseek(fid,0,'eof');
    fseek(fid,cPos,'bof');

    % ======================================
    %   constants and conditional variables
    % ======================================
        bit_0 = uint8(2^0);
        bit_5 = uint8(2^5);
        mdhStart = 1-byteMDH;
        
        u8_000 = zeros( 3, 1, 'uint8'); % for comparison with data_u8(1:3)

        % 20 fill bytes in VD (21:40)
        evIdx   = uint8(    21  + 20*isVD); % 1st byte of evalInfoMask
        dmaIdx  = uint8((29:32) + 20*isVD); % to correct DMA length using NCol and NCha
        if isVD
            dmaOff  = szScanHeader;
            dmaSkip = szChannelHeader;
        else
            dmaOff  = 0;
            dmaSkip = byteMDH;
        end
    % ======================================

    while true
        % Read mdh as binary (uint8) and evaluate as little as possible to know...
        %   ... where the next mdh is (ulDMALength / ushSamplesInScan & ushUsedChannels)
        %   ... whether it is only for sync (MDH_SYNCDATA)
        %   ... whether it is the last one (MDH_ACQEND)
        % evalMDH() contains the correct and readable code for all mdh entries.
        try
            % read everything and cut out the mdh
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
        catch exc
            warning( [mfilename() ':UnxpctdEOF'],  ...
                      [ '\nAn unexpected read error occurred at this byte offset: %d (%g GiB)\n'...
                        'Will stop reading now.\n'                                             ...
                        '=== MATLABs error message ================\n'                         ...
                        exc.message                                                            ...
                        '\n=== end of error =========================\n'                       ...
                       ], cPos, cPos/1024^3 )
            isEOF = true;
            break
        end

        bitMask = data_u8(evIdx);   % the initial 8 bit from evalInfoMask are enough

        if   isequal( data_u8(1:3), u8_000 )    ... % probably ulDMALength == 0
          || bitand(bitMask, bit_0)                % MDH_ACQEND

            % ok, look closer if really all *4* bytes are 0:
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );

            if ulDMALength == 0 || bitand(bitMask, bit_0)
                cPos = cPos + ulDMALength;
                % jump to next full 512 bytes
                if mod(cPos,512)
                    cPos = cPos + 512 - mod(cPos,512);
                end
                break;
            end
        end
        if bitand(bitMask, bit_5)  % MDH_SYNCDATA
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            cPos = cPos + ulDMALength;
            continue
        end

        % pehses: the pack bit indicates that multiple ADC are packed into one
        % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
        % since this code assumes one adc (x NCha) per DMA, we have to correct
        % the "DMA length"
        %     if mdh.ulPackBit
        % it seems that the packbit is not always set correctly
        NCol_NCha = double( typecast( data_u8(dmaIdx), 'uint16' ) );  % [ushSamplesInScan  ushUsedChannels]
        ulDMALength = dmaOff + (8*NCol_NCha(1) + dmaSkip) * NCol_NCha(2);

        n_acq = n_acq + 1;

        % grow arrays in batches
        if n_acq > szBlob
            filePos( end + allocSize ) = 0;
        end
        filePos( n_acq )  = cPos;

        cPos = cPos + ulDMALength;
    end % while true

    if isEOF
        n_acq = n_acq-1;    % ignore the last attempt
    end

    filePos( n_acq+1 ) = cPos;  % save pointer to the next scan
    filePos  = reshape( filePos(1:n_acq+1), 1, [] ); % row vector

end % of loop_mdh_read()

function change_twix_hdr(fid)
%% Reads in the distinct TWIX header parts, calls the search and replace
% function parse_deid.m, and overwrites the data directly in the file.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2016-08-25)
%
%       This code is heavily relying on TWIX parsing code bits from:
%       read_twix_hdr.m
%       (Dr. Philipp Ehses, University of Tuebingen)

        nbuffers = fread(fid, 1,'uint32');
        for b=1:nbuffers
            namesz  = 0;
            byte = 1;
            
            while byte~=0 % look for NULL-character
                byte   = fread(fid, 1, 'uint8');
                namesz = namesz+1;
            end
            
            fseek(fid,-namesz,'cof');
            bufname        = fread(fid,namesz,'char=>char').';
            bufname(end)   = []; % delete NULL character
            buflen         = fread(fid, 1,'uint32');
            buffer         = fread(fid, buflen, 'char=>char').';
            buff.(bufname) = buffer;
            
            % Search for patient info and replace with XXX happens now:
            newbuff.(bufname) = parse_deid(buff.(bufname));
            
            % Rewind by length of buffer and overwrite with new buffer
            fseek(fid,-buflen,'cof');
            fwrite(fid,newbuff.(bufname),'char');
        end
        
end % of change_twix_hdr()

function newbuffer = parse_deid(buffer)
%% This contains the regular expressions needed to filter out the
% information that requires de-identification. 
% 
% The code replaces the values assigned to the long/string parameter fields
%   PatientID
%   PatientBirthDay
%   PatientBirthDate
%   tPatientName
%   PatientsName
%   PatientName
%   lPatientSex
%   PatientSex
%   FrameOfReference (2017-02-02)
%   tFrameOfReference (2017-02-02)
%   FOR (2017-02-02)
%   ExamMemoryUID (2017-02-02)
%   tReferenceImage0 (2017-02-02)
%   tReferenceImage1 (2017-02-02)
%   tReferenceImage2 (2017-02-02)
%   DeviceSerialNumber (2017-02-06)
%   PatientLOID/PatientLoid (2018-09-14)
%   StudyLOID/StudyLoid (2018-09-14)
%   SeriesLOID (2018-09-14)
%   Study (2018-09-14)
%   Patient (2018-09-14)
%   InstitutionAddress (2018-09-14)
%   InstitutionName (2018-09-14)
%   StudyLOIDS4Split (2018-09-14)
%   MppsLOIDS4Split (2018-09-14)
%   SeriesLOIDS4Split (2018-09-14)
%   StudyLOIDS4SplitAndPause (2018-09-14)
%   MppsLOIDS4SplitAndPause (2018-09-14)
%   SeriesLOIDS4SplitAndPause (2018-09-14)
%   ParentLoid (2018-09-14)
%   ParentUid (2018-09-14)
%   DefaultSeriesLoid (2018-09-14)
%
% the double parameter fields with "precision" tag
%   flUsedPatientWeight
%   flPatientAge
%   PatientWeight
%   PatientAge
%
% the double parameter fields with "unit" tag
%   flPatientHeight
%
% the string parameter fields with "visible" tag
%   PatientName
%   PatientBirthDay
%   PatientSex
%   PatientID
%   PatientBirthDate
%   FrameOfReference (2017-02-02)
%   tFrameOfReference (2017-02-02)
%   FOR (2017-02-02)
%   ExamMemoryUID (2017-02-02)
%
% the double parameter fields with "visible" tag
%   PatientAge
%   PatientWeight
%
% the array parameter fields with "default" tag
%   ReferencedImageSequence (2017-02-02)
%   Studyloids4Split (2018-09-14)
%   MppsLoids4Split (2018-09-14)
%   DefaultSeriesLoids4Split (2018-09-14)
%
% the simple text field in the MeasYaps header part
%   tReferenceImage0 (2017-02-07)
%   tReferenceImage1 (2017-02-07)
%   tReferenceImage2 (2017-02-07)
%
% If your specific implementation requires more fields to be de-identified,
% please note the author.
%
% Author: Dr. Georg Oeltzschner
%         Johns Hopkins University, 08/23/2016
%
% History: 2016-08-23: First version.
%          2017-02-02: Added fields containing date stamps.
%          2017-02-07: Added fields containing device IDs.
%          2018-09-14: Added fields containing scan dates and unique IDs.

newbuffer1 = regexprep(buffer, '<Param(?:Bool|Long|String)\."(?:PatientID|PatientBirthDay|tPatientName|PatientName|PatientsName|PatientSex|lPatientSex|PatientBirthDate|FrameOfReference|tFrameOfReference|FOR|ExamMemoryUID|tReferenceImage0|tReferenceImage1|tReferenceImage2|DeviceSerialNumber|PatientLOID|PatientLoid|StudyLOID|StudyLoid|SeriesLOID|Study|Patient|InstitutionAddress|InstitutionName|ParentLoid|ParentUid|DefaultSeriesLoid)">\s*\n*\s*{\s*([^}]*)','${write_xxx($0,$1)}');
newbuffer2 = regexprep(newbuffer1, '<ParamDouble\."(?:flPatientAge|flUsedPatientWeight|PatientAge|PatientWeight)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','${write_xxx($0,$2)}');
newbuffer3 = regexprep(newbuffer2, '<ParamDouble\."(?:flPatientHeight)">\s*{\s*<Unit>.*"\[mm\]"\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','${write_xxx($0,$2)}');
newbuffer4 = regexprep(newbuffer3, '<Param(?:Bool|Long|String)\."(?:PatientName|PatientBirthDay|PatientSex|PatientID|PatientBirthDate|FrameOfReference|tFrameOfReference|FOR|ExamMemoryUID)">\s*{\s*(<Visible>\s*"\w*")?\s*([^}]*)','${write_xxx($0,$2)}');
newbuffer5 = regexprep(newbuffer4, '<ParamArray\."ReferencedImageSequence">\s*{\s*(<Default>\s*<ParamString."">)?\s*{\s*}\s*{\s*([^}]*)','${write_xxx($0,$2)}');
newbuffer6 = regexprep(newbuffer5, '<ParamArray\."ReferencedImageSequence">\s*{\s*(<Default>\s*<ParamString."">)?\s*{\s*}\s*{\s*([^}]*)}\s*{\s*([^}]*)}\s*{\s*([^}]*)','${write_xxx($0,$3)}');
newbuffer7 = regexprep(newbuffer6, '<ParamArray\."ReferencedImageSequence">\s*{\s*(<Default>\s*<ParamString."">)?\s*{\s*}\s*{\s*([^}]*)}\s*{\s*([^}]*)}\s*{\s*([^}]*)','${write_xxx($0,$4)}');
newbuffer8 = regexprep(newbuffer7,'"\d*\.\d*\.\d*\.\d*\.\d*\.\d*\.\d*\.\d*\.\d*\.\d*"','${write_xxx($0,$0)}');
newbuffer9 = regexprep(newbuffer8, '<Param(?:Bool|Long|String)\."(?:PatientID|PatientBirthDay|tPatientName|PatientName|PatientsName|PatientSex|lPatientSex|PatientBirthDate|FrameOfReference|tFrameOfReference|FOR|ExamMemoryUID|tReferenceImage0|tReferenceImage1|tReferenceImage2|DeviceSerialNumber|PatientLOID|StudyLOID|SeriesLOID|Study|Patient|InstitutionAddress|InstitutionName|StudyLOIDS4Split|MppsLOIDS4Split|SeriesLOIDS4Split|StudyLOIDS4SplitAndPause|MppsLOIDS4SplitAndPause|SeriesLOIDS4SplitAndPause)">\s*{\s*<MinSize>\s*[0-9]*\s*<MaxSize>\s*[0-9]*\s*([^}]*)','${write_xxx($0,$1)}');
newbuffer10 = regexprep(newbuffer9, '<ParamArray\."(?:Studyloids4Split|MppsLoids4Split|DefaultSeriesLoids4Split)">\s*{\s*<Default>\s*<ParamString."">\s*{\s*}\s*{\s*([^}]*)','${write_xxx($0,$1)}');
newbuffer11 = regexprep(newbuffer10, '<ParamArray\."(?:Studyloids4Split|MppsLoids4Split|DefaultSeriesLoids4Split)">\s*{\s*<MinSize>\s*[0-9]*\s*<MaxSize>\s*[0-9]*\s*<Default>\s*<ParamString."">\s*{\s*}\s*{\s*([^}]*)','${write_xxx($0,$1)}');

newbuffer = regexprep(newbuffer11, '<ParamDouble\."(?:PatientAge|PatientWeight)">\s*{\s*(<Visible>\s*"\w*")\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','${write_xxx($0,$3)}');

end % of parse_deid()

function [exitFunc, fnames] = CheckForOutput(nArgs, fnames)
%% Check if any de-identified files have already been output and ask user if
% they want to overwrite them
% Author:    
%       Dr. Mark Mikkelsen, Johns Hopkins University

exitFunc = 0;

if nArgs < 1
    
    if any(~cellfun('isempty', strfind(fnames, '_noID')))
        resp = input('\nDe-identified files found in the directory! Proceed and overwrite? [y/n]: ','s');
        if strcmpi(resp, 'y')
            disp('Overwriting...');
        elseif strcmpi(resp, 'n')
            disp('Exiting...');
            exitFunc = 1;
            return
        end
        
        ind1 = strfind(fnames, '_noID');
        ind2 = find(~cellfun('isempty', ind1));
        fnames(ind2) = []; %#ok<FNDSB>        
    end
    
else
    
    for ii = 1:length(fnames)
        fnames_noID{ii} = [fnames{ii}(1:end-5) '_noID' fnames{ii}(end-4:end)];
    end
    
    if any(cellfun(@exist, fnames_noID))
        resp = input('\nDe-identified files found in the directory! Proceed and overwrite? [y/n]: ','s');
        if strcmpi(resp, 'y')
            disp('Overwriting...');
        elseif strcmpi(resp, 'n')
            disp('Exiting...');
            exitFunc = 1;
            return
        end
        
    end
    
end

end % of CheckForOutput()