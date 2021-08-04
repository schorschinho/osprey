% op_preref.m
% USAGE:
% [refShift_ind_ini]=op_preref(in);
% 
% DESCRIPTION:
% We will use a frequency cross-correlation approach on the
% Choline and Creatine singlets to generate a robust initial
% frequency guess for the robust spectral registration. This is
% especially useful for data with heavy frequency drift. The
% transients are averaged into bins including 10% of the
% averages of the whole spectra and referenced afterwards. For
% these packages the same initial frequency guess is forwarded
% to op_robustSpecReg.
% 
% INPUTS:
% in	= input data in matlab structure format.  
%
% OUTPUTS:
% refShift_ind_ini   = Vector with initial frequency shifts
%
%
% AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2021-02-11)
%       hzoelln2@jhmi.edu
%
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-02-11: First version of the code.
function [refShift_ind_ini] = op_preref(raw, seq)

if nargin < 2
    seq = 'unedited';
end

% PREPARE DATA
% Based on the sequence type the raw data is split into the sub-spectra.
% Determine the polarity of the respective peak: if the absolute of the
% maximum minus the absolute of the minimum is positive, the polarity 
% of the respective peak is positive; if the absolute of the maximum 
% minus the absolute of the minimum is negative, the polarity is negative.
switch seq
    case 'unedited'
        temp_A = op_averaging(raw);
        raw_A_Cr    = op_freqrange(temp_A,2.8,3.2);
        polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
        temp_proc = raw;
        if polResidCr < 0        
            temp_proc = op_ampScale(temp_proc,-1);
        end
        
    case 'MEGA'
        % Get sub-spectra, depending on whether they are stored as such
        if raw.subspecs == 2
            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
        else
            raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
            raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
        end
        
        temp_A = op_averaging(raw_A);
        raw_A_Cr    = op_freqrange(temp_A,2.8,3.2);
        polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
        temp_rawA = raw_A;
        if polResidCr < 0        
            temp_rawA = op_ampScale(temp_rawA,-1);
        end

        temp_B = op_averaging(raw_B);
        raw_B_Cr    = op_freqrange(temp_B,2.8,3.2);
        polResidCr  = abs(max(real(raw_B_Cr.specs))) - abs(min(real(raw_B_Cr.specs)));
        temp_rawB = raw_B;
        if polResidCr < 0        
            temp_rawB = op_ampScale(temp_rawB,-1);
        end
        
        temp_proc   = op_addScans(temp_rawA,temp_rawB);
        
    case {'HERMES','HERCULES'}
        % Get sub-spectra, depending on whether they are stored as such
        if raw.subspecs == 4

            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
            raw_C   = op_takesubspec(raw,3);                    % Get third subspectrum
            raw_D   = op_takesubspec(raw,4);                    % Get fourth subspectrum

        else

            raw_A   = op_takeaverages(raw,1:4:raw.averages);    % Get first subspectrum
            raw_B   = op_takeaverages(raw,2:4:raw.averages);    % Get second subspectrum
            raw_C   = op_takeaverages(raw,3:4:raw.averages);    % Get third subspectrum
            raw_D   = op_takeaverages(raw,4:4:raw.averages);    % Get fourth subspectrum
        end

        temp_A = op_averaging(raw_A);
        raw_A_Cr    = op_freqrange(temp_A,2.8,3.2);

        polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
        temp_rawA = raw_A;
        if polResidCr < 0        
            temp_rawA = op_ampScale(temp_rawA,-1);
        end

        temp_B = op_averaging(raw_B);
        raw_B_Cr    = op_freqrange(temp_B,2.8,3.2);
        polResidCr  = abs(max(real(raw_B_Cr.specs))) - abs(min(real(raw_B_Cr.specs)));
        temp_rawB = raw_B;
        if polResidCr < 0        
            temp_rawB = op_ampScale(temp_rawB,-1);
        end

        temp_C = op_averaging(raw_C);
        raw_C_Cr    = op_freqrange(temp_C,2.8,3.2);
        polResidCr  = abs(max(real(raw_C_Cr.specs))) - abs(min(real(raw_C_Cr.specs)));
        temp_rawC = raw_C;
        if polResidCr < 0        
            temp_rawC = op_ampScale(temp_rawC,-1);
        end

        temp_D = op_averaging(raw_D);
        raw_D_Cr    = op_freqrange(temp_D,2.8,3.2);
        polResidCr  = abs(max(real(raw_D_Cr.specs))) - abs(min(real(raw_D_Cr.specs)));
        temp_rawD = raw_D;
        if polResidCr < 0        
            temp_rawD = op_ampScale(temp_rawD,-1);
        end
        temp_proc   = op_addScans(temp_rawA,temp_rawB);
        temp_proc   = op_addScans(temp_proc,temp_rawC);
        temp_proc   = op_addScans(temp_proc,temp_rawD);    
end


% CALCULATE FREQUENCY SHIFT
    % We will use a frequency cross-correlation approach on the
    % Choline and Creatine singlets to generate a robust initial
    % frequency guess for the robust spectral registration. This is
    % especially useful for data with heavy frequency drift. The
    % transients are averaged into bins including 10% of the
    % averages of the whole spectra and referenced afterwards. For
    % these packages the same initial frequency guess is forwarded
    % to op_robustSpecReg.
    
    temp_spec   = temp_proc;
    refShift_ind_ini = zeros(temp_proc.averages,1);
    for av = 1 : round(temp_proc.averages*0.1) :temp_proc.averages-(round(temp_proc.averages*0.1)-1)-mod(temp_proc.averages,round(temp_proc.averages*0.1)) % 10% packaging
        fids = temp_proc.fids(:,av:av+(round(temp_proc.averages*0.1)-1)); 
        specs = temp_proc.specs(:,av:av+(round(temp_proc.averages*0.1)-1));
        temp_spec.fids = mean(fids,2); % store average fid
        temp_spec.specs = mean(specs,2); % store average spectra
        [refShift, ~] = osp_CrChoReferencing(temp_spec); % determine frequency shift
        refShift_ind_ini(av : av+round(temp_proc.averages*0.1)-1) = refShift; %save initial frequency guess
    end
    if mod(temp_proc.averages,round(temp_proc.averages*0.1)) > 0 % remaining averages if data isn't a multiple of 10.
        fids = temp_proc.fids(:,end-(mod(temp_proc.averages,round(temp_proc.averages*0.1))-1):end);
        specs = temp_proc.specs(:,end-(mod(temp_proc.averages,round(temp_proc.averages*0.1))-1):end);
        temp_spec.fids = mean(fids,2); % store average fid
        temp_spec.specs = mean(specs,2); % store average spectra
        [refShift, ~] = osp_CrChoReferencing(temp_spec); % determine frequency shift
        refShift_ind_ini(end-(mod(temp_proc.averages,round(temp_proc.averages*0.1))-1) : temp_proc.averages) = refShift; %save initial frequency guess
    end
    
    
end