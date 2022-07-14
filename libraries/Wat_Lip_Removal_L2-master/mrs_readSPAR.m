function info = mrs_readSPAR( fileName )
% MRS_READSPAR reads Philips MRS header file (.SPAR) 
% info = mrs_readSPAR(fileName)
% ARGS :
% fileName = name of header file 
% RETURNS:
% info = header file information 
% EXAMPLE: 
% >> info = mrs_readSPAR('sub1.SPAR');
% >> disp(info)
% NOTE: 
% info.size: [LR,AP,FH]
% info.offcentre: [LR,AP,FH]
% info.angulation: [LR,AP,FH]
% info.FOV : [LR,AP,FH]
% info.dim : [LR,AP,FH]
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
% Copyright (c) 2013, University of Nottingham. All rights reserved.

	[~,~,ext]=fileparts(fileName);    
	if isempty(ext)==1
        fileName=[fileName,'.SPAR'];
	end
    
    fid=fopen(fileName,'r');
	header_info=textscan(fid,'%s','delimiter','\n');
	fclose(fid);
    
	no_lines=size(header_info{1});
    c=0;   
    info.CSI=0;
    for i=1:no_lines
        line=header_info{1}{i};
        s_ind = strfind(line,'samples');
        r_ind = strfind(line,'rows');
        tf_ind = strfind(line,'synthesizer_frequency');
        bw_ind = strfind(line, 'sample_frequency');
        te_ind=strfind(line, 'echo_time');
        tr_ind=strfind(line,'repetition_time');
        avg_ind = strfind(line, 'averages');
        
        ap_size_ind = strfind(line,'ap_size');
        lr_size_ind = strfind(line,'lr_size');
        cc_size_ind = strfind(line, 'cc_size');
        ap_off_centre_ind = strfind(line, 'ap_off_center');
        lr_off_centre_ind = strfind(line, 'lr_off_center');
        cc_off_centre_ind = strfind(line, 'cc_off_center');
        ap_angulation_ind = strfind(line, 'ap_angulation');
        lr_angulation_ind = strfind(line, 'lr_angulation');
        cc_angulation_ind = strfind(line, 'cc_angulation');
        phase_encode_ind = strfind(line, 'phase_encoding_enable');
        phase_encoding_ind=strfind(line, 'phase_encoding_fov');
        slice_thickness_ind=strfind(line, 'slice_thickness');
        
        
        dim1_ind = strfind(line,'dim2_pnts');
        dim2_ind = strfind(line,'dim3_pnts');
        dim3_ind = strfind(line,'nr_of_slices_for_multislice');
        
        
        if ~isempty(s_ind) % number of points in each spectrum
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.samples= str2double(str_temp{1}{3});
            
        elseif ~isempty(r_ind) % number of dynamics
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.dynamics = str2double(str_temp{1}{3});
            
        elseif ~isempty(tf_ind) % synthesizer frequency
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.transmit_frequency = str2double(str_temp{1}{3});
            
        elseif ~isempty(bw_ind) % spectral bandwidth
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.BW = str2double(str_temp{1}{3});
            
        elseif ~isempty(te_ind) % TE (echo time)
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.TE = str2double(str_temp{1}{3});
            
        elseif ~isempty(tr_ind) % TE (echo time)
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.TR = str2double(str_temp{1}{3});
            
        elseif ~isempty(avg_ind) % number of spectra averaged over
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.averages = str2double(str_temp{1}{3});
            
        elseif ~isempty(ap_size_ind) % Anterior-Posterior size of VOI
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.size(2)= str2double(str_temp{1}{3});
            
        elseif ~isempty(lr_size_ind) % Left-Right size of VOI
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.size(1) = str2double(str_temp{1}{3});
            
        elseif ~isempty(cc_size_ind) % Feet-Head size of VOI
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.size(3) = str2double(str_temp{1}{3});
            
        elseif ~isempty(ap_off_centre_ind) % off centre (AP) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.offcentre(2) = str2double(str_temp{1}{3});
            
        elseif ~isempty(lr_off_centre_ind) % off centre (LR) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.offcentre(1) = str2double(str_temp{1}{3});
            
        elseif ~isempty(cc_off_centre_ind) % off centre (FH) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.offcentre(3) = str2double(str_temp{1}{3});
            
        elseif ~isempty(ap_angulation_ind) % angulation (AP) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.angulation(2) = str2double(str_temp{1}{3});
            
        elseif ~isempty(lr_angulation_ind) % angulation (LR) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.angulation(1) = str2double(str_temp{1}{3});
            
        elseif ~isempty(cc_angulation_ind) % angulation (FH) of VOI or FOV
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.angulation(3) = str2double(str_temp{1}{3});
            
        elseif ~isempty(phase_encode_ind) % phase encoded?
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            if length(str_temp{1}{3}) == 5;
                info.CSI=1;
            end
            
        elseif ~isempty(phase_encoding_ind) % FOV, x, LR
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.FOV(1) = str2double(str_temp{1}{3});
            
        elseif ~isempty(slice_thickness_ind) % slice thickness, z, FH         
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            if c==0;
                info.FOV(3) = str2double(str_temp{1}{3});
                c=1;
            end
            
        elseif ~isempty(dim1_ind) % dimension 1, x, RL 
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.dim(1) = str2double(str_temp{1}{3});
            
        elseif ~isempty(dim2_ind) % dimension 2, y, AP
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.dim(2) = str2double(str_temp{1}{3});
            
        elseif ~isempty(dim3_ind) % dimsion 3, z, FH
            str_temp = textscan(line, '%s', 'delimiter', ' ');
            info.dim(3) = str2double(str_temp{1}{3});
        end
    end
    if info.CSI==1
       info.FOV(2)=info.FOV(1)/info.dim(1)*info.dim(2); % y, AP
    end
end

