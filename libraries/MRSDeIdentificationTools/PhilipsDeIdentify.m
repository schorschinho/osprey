function PhilipsDeIdentify(fnames)
% PhilipsDeIdentify
%   Reads Philips SPAR files and removes participant information. Newly
%   de-identified SPAR files are then output, with filenames appended with
%   '_noID'. The original files are not overwritten.
%
%   PhilipsDeIdentify will also create copies of the SDAT files associated
%   with the SPAR files, with filenames appended with '_noID'. This step is
%   required for Gannet functionality. (It is assumed that each SDAT/SPAR
%   file pair has the same name.)
%
%   If PhilipsDeIdentify has already been run, and de-identified files are
%   found within the directory, the user will be asked whether these files
%   should be overwritten.
%
%   NOTE: The user must make sure that filenames themselves do not contain
%   information that can personally identify participants. This function
%   will only de-identify the content of the SPAR file.
%
%   Usage:
%       PhilipsDeIdentify, by itself, de-identifies all SPAR files found
%       within the current directory.
%
%       PhilipsDeIdentify(fnames) de-identifies SPAR files listed in the
%       cell array fnames.
%
%   Example:
%       c = {'S01_gaba_7_2_raw_act.SPAR', 'S01_gaba_7_2_raw_ref.SPAR'};
%       PhilipsDeIdentify(c);

%   Author: Mark Mikkelsen (Johns Hopkins University, 2018)
%
%   Version history:
%       2016-08-03: + Function created
%       2016-08-04: + All files found in current directory de-identified
%                     if nargin < 1
%                   + Copies of SDAT files created; filenames appended with
%                     '_noID' (req. for Gannet)
%                   + CheckForOutput added
%       2017-01-19: + Fix for case-sensitivity of extensions
%       2018-09-13: + Remove scan date and time
%       2018-09-25: + Minor bug fix


if nargin < 1 % De-identify all SPAR files in current directory
    
    flist = dir('*.spar');
    if isempty(flist)
        flist = dir('*.SPAR');
    end
    flist = flist(cellfun(@isempty, strfind({flist.name}, '._'))); %#ok<*STRCLFH>
    for ii = 1:length(flist)
        fnames(ii) = cellstr(flist(ii).name);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
else % De-identify SPAR files user has listed in fnames
    
    % Check if filenames include a .SPAR/.spar extension
    for ii = 1:length(fnames)
        ext = fnames{ii}(end-4:end);
        assert(strcmpi(ext, '.spar'), ...
            ['The filename ' fnames{ii} ' does not include a .SPAR/.spar extension.']);
    end
    
    % Check if files can be found
    for ii = 1:length(fnames)
        assert(any(exist(fnames{ii}, 'file')), ...
            ['The file ' fnames{ii} ' cannot be found.' ...
            ' Check spelling of filenames (SDAT/SPAR files must include an extension in their filename).' ...
            ' Also check that you are in the right directory.']);
    end
    
    nArgs = nargin;
    [exitFunc, fnames] = CheckForOutput(nArgs, fnames);
    if exitFunc
        return
    end
    
end

% Read SPAR files and remove participant information; save new SPAR files
for ii = 1:length(fnames)
    
    spar_fid = fopen(fnames{ii}, 'r');
    spar_fid_noID = fopen([fnames{ii}(1:end-5) '_noID' fnames{ii}(end-4:end)], 'w');
    
    tline = fgetl(spar_fid);
    while ischar(tline)
        if any(strfind(tline, 'examination_name'))
            tline = 'examination_name : ';
        elseif any(strfind(tline, 'patient_name'))
            tline = 'patient_name : ';
        elseif any(strfind(tline, 'patient_birth_date'))
            tline = 'patient_birth_date : ';
        elseif any(strfind(tline, 'scan_date'))
            tline = 'scan_date : ';
        end
        fprintf(spar_fid_noID, '%s\n', tline);
        tline = fgetl(spar_fid);
    end
    
    fclose(spar_fid);
    fclose(spar_fid_noID);
    
end

% Create duplicate SDAT files
for ii = 1:length(fnames)
    
    if all(isstrprop(fnames{ii}(end-3:end), 'lower'))
        fnames_sdat = [fnames{ii}(1:end-5) '.sdat']; %#ok<*AGROW>
    elseif all(isstrprop(fnames{ii}(end-3:end), 'upper'))
        fnames_sdat = [fnames{ii}(1:end-5) '.SDAT'];
    else
        fnames_sdat = [fnames{ii}(1:end-5) '.sdat'];
    end
    
    assert(any(exist(fnames_sdat, 'file')), ...
        ['The file ' fnames_sdat ' cannot be found.' ...
        ' SDAT/SPAR file pairs must have the same name (and include an extension).' ...
        ' Also check that all files are in the same directory.']);
    copyfile(fnames_sdat, [fnames_sdat(1:end-5) '_noID' fnames_sdat(end-4:end)]);
    
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



