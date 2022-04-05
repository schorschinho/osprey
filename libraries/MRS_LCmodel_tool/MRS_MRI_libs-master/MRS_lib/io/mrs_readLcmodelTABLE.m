function info = mrs_readLcmodelTABLE( fileName )
% MRS_READLCMODELTABLE reads the important metabolites concentration from the 
% LCModel output file (.table)
%
% info = mrs_readLcmodelTABLE( fileName )
%
% ARGS :
% fileName = name of LCModel ouput .table file 
%
% RETURNS:
% info = information about quantified metabolite concentration, standard
% deviation, and relative concentration with respective to Pcr+Cr.
% 
% EXAMPLE: 
% >> info = mrs_readLcmodelTABLE('sub6_sl1_1-1.table');
% >> info.name % name of the metabolites quantified 
% >> info.concentration % absolute concentration of the metabolites quantified
% >> info.SDpct % SD% of the metabolites quantified
% >> info.relative_conc % relative concentration of the metabolites quantified
%
%
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
%
% Copyright (c) 2013, University of Nottingham. All rights reserved.
    
	[~,~,ext]=fileparts(fileName);    
    if isempty(ext)==1
        fileName=[fileName,'.table'];
    end
    
    fid=fopen(fileName,'r');
    header_info=textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    
    no_lines=size(header_info{1});
    metabolites_num=0;
    for i=1:no_lines
        line=header_info{1}{i};    
        
        if ~isempty(strfind(line, '$$CONC'))
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            metabolites_num=str2double(str_temp{1}{2})-1;            
        end
        
        if ~isempty(strfind(line, 'Conc.'))
            break;
        end
    end
    n=0;
    for m=i+(1:metabolites_num)
        n=n+1;
        line=header_info{1}{m};
        str_temp = textscan(line, '%s', 'delimiter', ' ');
        r=0;
        temp={};
        for c=1:length(str_temp{1})
            if ~isempty(str_temp{1}{c})
                r=r+1;
                temp{r}=str_temp{1}{c};
            end
        end
        if length(temp)==4 
            info.name(n)=temp(4);
            info.concentration(n)=str2double(temp{1});
            sd=temp{2};
            info.SDpct(n)=str2double(sd(1:end-1));
            info.relative_conc(n)=str2double(temp{3});            
        else
            if ~isempty(strfind(line, '+'))
                l=textscan(temp{3}, '%s', 'delimiter', '+');
            elseif ~isempty(strfind(line, '-'))
                l=textscan(temp{3}, '%s', 'delimiter', '-');
            end
            info.name(n)=l{1}(2);
            info.concentration(n)=str2double(temp{1});
            sd=temp{2};
            info.SDpct(n)=str2double(sd(1:end-1));
            info.relative_conc(n)=str2double(l{1}{1}); 
        end
    end
    
    % Extract additional modeling parameters from $$MISC bit
    % Find the line containing the $$MISC string
    indexMisc = find(~cellfun(@isempty, strfind(header_info{1}, '$$MISC')));
    
    % Get the FWHM/SNR line
    line = header_info{1}{indexMisc+1};
    str_temp = textscan(line, '%s %s %f %s %s %s %f');
    info.fwhm = str_temp{3};
    info.snr = str_temp{7};
    
    % Get the data shift line
    line = header_info{1}{indexMisc+2};
    str_temp = textscan(line, '%s', 'delimiter', '=');
    str_ppm = textscan(str_temp{1}{2}, '%f %s');
    info.refShift = str_ppm{1};
    
    % Get the phase line
    line = header_info{1}{indexMisc+3};
    str_temp = textscan(line, '%s %f %s %f %s');
    info.ph0 = str_temp{2};
    info.ph1 = str_temp{4};
        

    
 end
      

