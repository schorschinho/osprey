function GEDeIdentify(fnames)
% GEDeIdentify
%   Reads GE P-files and removes participant/exam information from the
%   header. Newly de-identified P-files are then output, with filenames
%   appended with '_noID'. The original files are not overwritten.
%
%   If GEDeIdentify has already been run, and de-identified files are found
%   within the directory, the user will be asked whether these files should
%   be overwritten.
%
%   NOTE: The user must make sure that filenames themselves do not contain
%   information that can personally identify participants. This function
%   will only de-identify the content of the P-file header.
%
%   Usage:
%       GEDeIdentify, by itself, de-identifies all P-files found within the
%       current directory.
%
%       GEDeIdentify(fnames) de-identifies P-files listed in the cell array
%       fnames.
%
%   Example:
%       c = {'P15221.7', 'P17981.7'};
%       GEDeIdentify(c);

%   Author: Mark Mikkelsen (Johns Hopkins University, 2018)
%
%   Acknowledgements:
%       This code is heavily based on plotp.m and pfile_anon.c written by
%       Fred Frigo (Marquette University). The author is also grateful to
%       Ralph Noeske (GE Healthcare) for providing details regarding the
%       P-file header.
%
%   Version history:
%       2016-08-09: + Function created
%       2018-09-05: + Smarter parsing of P-file header
%       2018-09-25: + Minor bug fix
%       2020-10-23: + Added support for rdbm_rev_num 27.x

if nargin < 1 % De-identify all P-files in current directory
    
    flist = dir('*.7');
    flist = flist(cellfun(@isempty, strfind({flist.name}, '._'))); %#ok<*STRCLFH>
    for ii = 1:length(flist)
        fnames(ii) = cellstr(flist(ii).name);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
else % De-identify P-files user has listed in fnames
    
    % Check if filenames include a .7 extension
    for ii = 1:length(fnames)
        ext = fnames{ii}(end-1:end);
        assert(strcmpi(ext, '.7'), ['The filename ' fnames{ii} ' does not include a .7 extension.']);
    end
    
    % Check if files can be found
    for ii = 1:length(fnames)
        assert(any(exist(fnames{ii}, 'file')), ...
            ['The file ' fnames{ii} ' cannot be found.' ...
            ' Check spelling of filenames (P-files must include an extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
end

% Read P-files and remove participant/exam information; save new P-files
for ii = 1:length(fnames)
    
    pfile_fid = fopen(fnames{ii}, 'r', 'ieee-be');
    pfile_fid_noID = fopen([fnames{ii}(1:end-2) '_noID' fnames{ii}(end-1:end)], 'w', 'ieee-be');
    
    % Check if P-file can be found
    if pfile_fid == -1
        fclose(pfile_fid);
        fclose(pfile_fid_noID);
        error(['The file ' fnames{ii} ' cannot be found.' ...
            ' Check spelling of filenames (P-files must include an extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end
    
    % Determine P-file version
    [pfile_fid, pfile_fid_noID, hdr] = VersionCheck(fnames{ii}, pfile_fid, pfile_fid_noID);
    
    % Determine how big the P-file is
    file_info = dir(fnames{ii});
    pfile_sz = file_info.bytes;
    
    % Read the entire P-file
    frewind(pfile_fid);
    pfile = fread(pfile_fid, pfile_sz, 'integer*2');
    
    % Spectro prescan P-files
    if hdr.n_points == 1 && hdr.n_rows == 1
        hdr.n_points = 2048;
    end
    
    % Compute size (in bytes) of raw data
    data_elements = hdr.n_points * 2;
    total_frames = hdr.n_rows * hdr.n_echoes;
    data_elements = data_elements * total_frames * hdr.n_rcvrs;
    
    % Read the raw data (using the appropriate precision)
    fseek(pfile_fid, hdr.hdr_sz, 'bof');
    if hdr.point_sz == 2
        pfile_raw_data = fread(pfile_fid, data_elements, 'integer*2');
    else
        pfile_raw_data = fread(pfile_fid, data_elements, 'integer*4');
    end
    
    % Write all P-file bytes into new P-file
    frewind(pfile_fid_noID);
    fwrite(pfile_fid_noID, pfile, 'integer*2');
    
    % Overwrite EXAMDATATYPE structure that contains participant/exam
    % information with zeros in new P-file
    exam_data_type_sz = hdr.series_offset - hdr.exam_offset;
    fseek(pfile_fid_noID, hdr.exam_offset, 'bof');
    fwrite(pfile_fid_noID, zeros(exam_data_type_sz,1), 'integer*2');
    
    % Write raw data to new P-file (using the appropriate precision)
    fseek(pfile_fid_noID, hdr.hdr_sz, 'bof');
    if hdr.point_sz == 2
        fwrite(pfile_fid_noID, pfile_raw_data, 'integer*2');
    else
        fwrite(pfile_fid_noID, pfile_raw_data, 'integer*4');
    end
    
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    
end


function [exitFunc, fnames] = CheckForOutput(nArgs, fnames)
% Check if any de-identified files have already been output and ask user if
% they want to overwrite them

exitFunc = 0;

if nArgs < 1
    
    if any(~cellfun('isempty', strfind(fnames, '_noID'))) %#ok<*STRCL1>
        resp = input('\nDe-identified files found in the directory! Proceed and overwrite? [y/n]: ','s');
        if strcmpi(resp, 'y')
            disp('Overwriting...');
        elseif strcmpi(resp, 'n')
            disp('Exiting...');
            exitFunc = 1;
            return
        end
        
        ind1 = strfind(fnames, '_noID');
        ind2 = ~cellfun('isempty', ind1);
        fnames(ind2) = [];
    end
    
else
    
    for ii = 1:length(fnames)
        fnames_noID{ii} = [fnames{ii}(1:end-2) '_noID' fnames{ii}(end-1:end)]; %#ok<AGROW>
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


function [pfile_fid, pfile_fid_noID, hdr] = VersionCheck(fnames, pfile_fid, pfile_fid_noID)
% Determine P-file version

frewind(pfile_fid);
rdbm_rev_num = fread(pfile_fid, 1, 'real*4');

if rdbm_rev_num == 7.0 % Signa 8.X (LX)
    hdr.hdr_sz = 39984;
    hdr.exam_offset = hdr.hdr_sz - 1040 - 1028 - 1044;
    hdr.series_offset = hdr.hdr_sz - 1028 - 1044;
elseif rdbm_rev_num == 8.0 % Cardiac / MGD
    %hdr.hdr_sz = 60464;
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    error('GEDeIdentify does not yet support P-file header revision: %g', rdbm_rev_num);
elseif rdbm_rev_num > 5.0 && rdbm_rev_num < 6.0 % Signa 5.5
    %hdr.hdr_sz = 39940;
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    error('GEDeIdentify does not yet support P-file header revision: %g', rdbm_rev_num);
else
    % In 11.0 and later, the header and data are stored as little-endian
    fclose(pfile_fid);
    fclose(pfile_fid_noID);
    pfile_fid = fopen(fnames, 'r', 'ieee-le');
    pfile_fid_noID = fopen([fnames(1:end-2) '_noID' fnames(end-1:end)], 'w', 'ieee-le');
    
    frewind(pfile_fid);
    rdbm_rev_num = fread(pfile_fid, 1, 'real*4');
    
    if rdbm_rev_num == 9.0  % Excite2 11.0 product release
        hdr.hdr_sz = 61464;
        hdr.exam_offset = hdr.hdr_sz - 1040 - 1536 - 1536;
        hdr.series_offset = hdr.hdr_sz - 1536 - 1536;
    elseif rdbm_rev_num == 11.0  % 12.0 product release
        fseek(pfile_fid, 1468, 'bof');
        hdr.hdr_sz = fread(pfile_fid, 1, 'integer*4');
        fseek(pfile_fid, 1496, 'bof');
        hdr.exam_offset = fread(pfile_fid, 1, 'integer*4');
        fseek(pfile_fid, 1500, 'bof');
        hdr.series_offset = fread(pfile_fid, 1, 'integer*4');
    elseif rdbm_rev_num > 11.0
        chkRev = {'14.3','16','20.006','20.007','24','26.002','27','27.001'};
        if ~any(strcmp(num2str(rdbm_rev_num), chkRev))
            fclose(pfile_fid);
            fclose(pfile_fid_noID);
            error('GEDeIdentify does not yet support P-file header revision: %g', rdbm_rev_num);
        end
        
        switch num2str(rdbm_rev_num)
            case {'14.3','16','20.006','20.007','24'}
                rdb_hdr_off_image     = 377;
                rdb_hdr_off_data      = 368;
                rdb_hdr_off_exam      = 375;
                rdb_hdr_off_series    = 376;
                rdb_hdr_ps_mps_freq   = 107;
                rdb_hdr_nechoes       = 36;
                rdb_hdr_point_size    = 42;
                rdb_hdr_da_xres       = 52;
                rdb_hdr_da_yres       = 53;
                rdb_hdr_dab_start_rcv = 101;
                rdb_hdr_dab_stop_rcv  = 102;
            case {'26.002','27','27.001'}
                rdb_hdr_off_image     = 11;
                rdb_hdr_off_data      = 2;
                rdb_hdr_off_exam      = 9;
                rdb_hdr_off_series    = 10;
                rdb_hdr_ps_mps_freq   = 123;
                rdb_hdr_nechoes       = 74;
                rdb_hdr_point_size    = 80;
                rdb_hdr_da_xres       = 90;
                rdb_hdr_da_yres       = 91;
                rdb_hdr_dab_start_rcv = 133;
                rdb_hdr_dab_stop_rcv  = 134;
        end
        
        frewind(pfile_fid);
        i_hdr_value = fread(pfile_fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');
        frewind(pfile_fid);
        hdr_value = fread(pfile_fid, rdb_hdr_dab_stop_rcv, 'integer*2');
        
        hdr.hdr_sz        = i_hdr_value(rdb_hdr_off_data);
        hdr.exam_offset   = i_hdr_value(rdb_hdr_off_exam);
        hdr.series_offset = i_hdr_value(rdb_hdr_off_series);
        
        hdr.n_echoes = hdr_value(rdb_hdr_nechoes);
        hdr.point_sz = hdr_value(rdb_hdr_point_size);
        hdr.n_points = hdr_value(rdb_hdr_da_xres);
        hdr.n_rows   = hdr_value(rdb_hdr_da_yres);
        start_recv   = hdr_value(rdb_hdr_dab_start_rcv);
        stop_recv    = hdr_value(rdb_hdr_dab_stop_rcv);
        hdr.n_rcvrs  = (stop_recv - start_recv) + 1;
    end
end



