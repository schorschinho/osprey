function [info ] = mrs_readLcmodelPRINT( fileName )
%MRS_READLCMODELPRINT reads the LCmodel output file .print, which contains the
%all infos about the modelling process.
%
% [info ] = mrs_readLcmodelPRINT( fileName )
%
% INPUT :
% fileName = name of LCModel output .print file 
%
% OUTPUT:
% info = stores some interesting infos
% 
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2019-02-19)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-11-11: First version of the code.

    [~,~,ext]=fileparts(fileName);  
    
    if isempty(ext)==1
        fileName=[fileName,'.print'];
    end
    
    fid=fopen(fileName,'r');
	f=textscan(fid,'%s','delimiter','\n');
	fclose(fid);
    
	no_lines=size(f{1});
    info=[];
    for i=1:no_lines
        line=f{1}{i};   
        inital_ind = strfind(line, 'Starting values for final analysis');       
        if ~isempty(inital_ind)
                line=f{1}{i+2};
                str_temp = textscan(line, '%s', 'delimiter', ' ');               
                info.iniph0= str2double(str_temp{1}{end-1});
                line=f{1}{i+3};
                str_temp = textscan(line, '%s', 'delimiter', ' ');               
                info.iniph1= str2double(str_temp{1}{end-1});
                line=f{1}{i+4};
                str_temp = textscan(line, '%s', 'delimiter', ' ');               
                info.iniFWHM= str2double(str_temp{1}{3});
        end        
    end
    
    for i=1:no_lines
        line=f{1}{i};
        inital_ind = strfind(line, 'Area of unsuppressed water peak');
        if ~isempty(inital_ind)
            str_temp = textscan(line, '%s', 'delimiter', '=');
            info.h2oarea = str2double(str_temp{1}{2});  
        end
    end
    

end

