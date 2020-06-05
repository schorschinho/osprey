function [ spectra, spectra_metabolites, x_ppm, info ] = mrs_readLcmodelCOORD( fileName )
%MRS_READLCMODELCOORD reads the LCmodel output file .coord, which contains the
%coordinates of all curves on the one-page output.
%
% [ spectra, spectra_metabolites, x_ppm, info ] = mrs_readLcmodelCOORD(fileName)
%
% ARGS :
% fileName = name of LCModel output .coord file 
%
% RETURNS:
% spectra = matrix contains data spectrum, fitted spectrum, baseline [only the real values]    
% spectra_metabolites = matrix contains fitted spectrum for each metabolite with baseline [only the real values]    
% x_ppm = frequency axis in ppm
% info = info.n (points in each spectrum); info.metabolites (name list of 
%        all fitted metabolites)
% 
% EXAMPLE: 
% >>[spectra, spectra_metabolites, x_ppm, info]=mrs_readLcmodelCOORD('sub1_MRS.coord');
% >>figure; plot(x_ppm,spectra); set(gca,'XDir','reverse');
% >>legend('data','fitted spectrum','baseline');
% >>disp(info.metabolites); %This gives the order of the metabolite spectra
% in spectra_metabolites. 
% >> baseline=spectra(:,3);
% >> figure; plot(x_ppm,spectra_metabolites-baseline); set(gca,'XDir','reverse');
% >> legend(info.metabolites,'Location','EastOutside');
%
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
%
% Copyright (c) 2013, University of Nottingham. All rights reserved.

    [~,~,ext]=fileparts(fileName);  
    
    if isempty(ext)==1
        fileName=[fileName,'.coord'];
    end
    
    fid=fopen(fileName,'r');
	f=textscan(fid,'%s','delimiter','\n');
	fclose(fid);
    
	no_lines=size(f{1});
    m=0; 
    l=[];
    for i=1:no_lines
        line=f{1}{i};   
        ppm_ind = strfind(line, 'ppm-axis');
        conc_ind = strfind(line, 'Conc. =');
        ph_ind = strfind(line, 'Ph:');        
        if ~isempty(ppm_ind)
                str_temp = textscan(line, '%s', 'delimiter', ' ');               
                info.n= str2double(str_temp{1}{1});
                l=[l,i];
        end
        
        if ~isempty(conc_ind)
                m=m+1;
                str_temp = textscan(line, '%s', 'delimiter', ' ');               
                info.metabolites{m}=str_temp{1}{1}; 
                l=[l,i];
        end 
        
        if ~isempty(ph_ind)
                str_temp = textscan(line, '%s', 'delimiter', ' '); 
                if ~isempty(str_temp{1}{4}) && ~isnan(str2double(str_temp{1}{4}))
                    info.ph0= str2double(str_temp{1}{4});
                else
                    if ~isempty(str_temp{1}{3}) && ~isnan(str2double(str_temp{1}{3}))
                    info.ph0= str2double(str_temp{1}{3});
                    else
                    info.ph0= str2double(str_temp{1}{2});  
                    end
                end
                if ~isempty(str_temp{1}{12}) && ~isnan(str2double(str_temp{1}{12}))
                    info.ph1= str2double(str_temp{1}{12});
                else
                    if ~isempty(str_temp{1}{11}) && ~isnan(str2double(str_temp{1}{11}))
                    info.ph1= str2double(str_temp{1}{11});
                    else
                    if ~isempty(str_temp{1}{10}) && ~isnan(str2double(str_temp{1}{10}))
                        
                        info.ph1= str2double(str_temp{1}{10}); 
                    else
                        info.ph1= str2double(str_temp{1}{9}); 
                    end
                    end
                end
        end
    end
   
    if ~isempty(l)        
        % x axis (ppm)
        fid = fopen(fileName);
        C = textscan(fid, '%f', info.n, 'headerlines',l(1));
        x_ppm = C{1,1};
        fclose(fid);
        
        if length(l)>1 
            % data spectrum, fitted spetrum, baseline
            slist=l(1)+(ceil(info.n/10)+1):(ceil(info.n/10)+1):l(2)-(ceil(info.n/10)+1);
            for i=1:length(slist)
                fid = fopen(fileName);
                C = textscan(fid, '%f', info.n, 'headerlines', slist(i));
                spectra(:,i) = C{1,1};
                fclose(fid);
            end
            % individual curves for fitted metabolite spectrum
            for i=2:length(l)
                fid = fopen(fileName);
                C = textscan(fid, '%f', info.n, 'headerlines', l(i));
                spectra_metabolites(:,i-1) = C{1,1};
                fclose(fid);           
            end
        end
    end

end

