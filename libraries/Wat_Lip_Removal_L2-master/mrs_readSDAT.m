function data = mrs_readSDAT( fileName )
% MRS_READSDAT reads Philips MRS data file (.SDAT) 
% data = mrs_readSDAT(fileName)
% ARGS :
% fileName = name of data file 
% RETURNS:
% data = averaged FIDs  
% EXAMPLE: 
% >> fid = mrs_readSDAT('sub1.SDAT');
% >> size(fid)
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
% Copyright (c) 2013, University of Nottingham. All rights reserved.
    
    [~,file_name,ext]=fileparts(fileName);  
    info=mrs_readSPAR([file_name,'.SPAR']);
    
    if isempty(ext)==1
        fileName=[fileName,'.SDAT'];
    end
    
    % get date length
    fid = fopen(fileName,'r','ieee-le');
	data_size=length(fread(fid));
	fclose(fid);
            
	% read in data
	fid = fopen(fileName,'r','ieee-le');
	data=freadVAXG(fid,data_size,'float32');
	fclose(fid);
	
    data=reshape(data,2,[]);    
    data=data(1,:)+1i*data(2,:);
   
    if info.CSI==0
        data=reshape(data,info.samples,[]);     
    else
        data=squeeze(reshape(data,info.samples,info.dim(1),info.dim(2),[]));
    end
    
end

