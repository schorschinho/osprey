function [MRSCont] = osp_combineCoils(MRSCont, kk, ll, ref_ll, w_ll)
%% [MRSCont] = osp_combineCoils(MRSCont)
%   This function performs a the receiver coil combination of multi-array
%   data. All coil-combination procedures are performed using the ratio of
%   the maximum signal in each receiver to the square of the noise as the
%   weighting factor (Hall et al., NeuroImage 86:35-42 (2014)).
%
%   If the MRSCont structure contains a reference scan (i.e. data
%   acquired with the same TE and sequence as the metabolite data), the
%   metabolite and reference data are combined based on this scan. If there
%   is no reference scan, the metabolite data is combined based on its own
%   coil sensitivities.
%
%   If MRSCont contains a (short-TE) water scan, it is combined separately
%   using coil sensitivities derived from its own signals.
%
%   USAGE:
%       [MRSCont] = osp_combineCoils(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-20)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-20: First version of the code.

% For SPECIAL acquisitions, some of the sub-spectra need to be combined
% prior to determining the CC coefficients. We'll set a flag here.
isSpecial = strcmpi(MRSCont.raw_uncomb{kk}.seq, 'special');
MRSCont.flags.isSPECIAL = isSpecial;

%% Calculate coil combination weights
if nargin<5
    
    % Loop over all datasets
    for kk = 1:MRSCont.nDatasets(1)
        for ll = 1: 1:MRSCont.nDatasets(2)
            
            % Check if reference scans exist, if so, get CC coefficients from there
            if MRSCont.flags.hasRef
                
                if isSpecial
                    % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
                    cweights          = op_getcoilcombos(op_combinesubspecs(MRSCont.raw_ref_uncomb{ll,kk}, 'diff'), 1, 'h');
                else
                    cweights          = op_getcoilcombos(MRSCont.raw_ref_uncomb{ll,kk},1,'h');
                end
                
                raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{ll,kk},1,'h',cweights);
                raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{ll,kk},1,'h',cweights);
                
                if MRSCont.flags.isUnEdited
                    raw_comb.flags.isUnEdited = 1;
                    raw_ref_comb.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA
                    raw_comb.flags.isMEGA = 1;
                    raw_ref_comb.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES
                    raw_comb.flags.isHERMES = 1;
                    raw_ref_comb.flags.isHERMES = 1;
                elseif MRSCont.flags.isHERCULES
                    raw_comb.flags.isHERCULES = 1;
                    raw_ref_comb.flags.isHERCULES = 1;
                end
                
                MRSCont.raw{ll,kk}     = raw_comb;
                MRSCont.raw_ref{ll,kk} = raw_ref_comb;
                if MRSCont.raw_ref{ll,kk}.subspecs > 1 && ~isSpecial
                    if MRSCont.flags.isMEGA
                        raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ll,kk},1);
                        raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ll,kk},2);
                        MRSCont.raw_ref{ll,kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                    else
                        raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ll,kk},1);
                        raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ll,kk},2);
                        raw_ref_C               = op_takesubspec(MRSCont.raw_ref{ll,kk},3);
                        raw_ref_D               = op_takesubspec(MRSCont.raw_ref{ll,kk},4);
                        MRSCont.raw_ref{ll,kk} = op_concatAverages(raw_ref_A,raw_ref_B,raw_ref_C,raw_ref_D);
                    end
                end
                
            else
                
                % if not, use the metabolite scan itself
                if isSpecial
                    % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
                    cweights          = op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(MRSCont.raw_uncomb{ll,kk}), 'diff'), 0, 0.01, 2);
                else
                    cweights          = op_getcoilcombos(MRSCont.raw_uncomb{ll,kk}, 1, 'h');
                end
                
                raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{ll,kk},1,'h',cweights);
                if MRSCont.flags.isUnEdited
                    raw_comb.flags.isUnEdited = 1;
                elseif MRSCont.flags.isMEGA
                    raw_comb.flags.isMEGA = 1;
                elseif MRSCont.flags.isHERMES
                    raw_comb.flags.isHERMES = 1;
                elseif MRSCont.flags.isHERCULES
                    raw_comb.flags.isHERCULES = 1;
                end
                MRSCont.raw{ll,kk}     = raw_comb;
            end
            
            % Now do the same for the (short-TE) water signal
            if MRSCont.flags.hasWater
                if isSpecial
                    % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
                    cweights_w          = op_getcoilcombos(op_combinesubspecs(MRSCont.raw_w_uncomb{ll,kk}, 'diff'), 1, 'h');
                else
                    cweights_w          = op_getcoilcombos(MRSCont.raw_w_uncomb{ll,kk}, 1, 'h');
                end
                raw_w_comb          = op_addrcvrs(MRSCont.raw_w_uncomb{ll,kk},1,'h',cweights_w);
                raw_w_comb.flags.isUnEdited = 1;
                MRSCont.raw_w{ll,kk}   = raw_w_comb;
            end
        end
        
    end
    
else
    
    % For SPECIAL acquisitions, some of the sub-spectra need to be combined
    % prior to determining the CC coefficients. We'll set a flag here.
    isSpecial = strcmpi(MRSCont.raw_uncomb{kk}.seq, 'special');
    
    % Check if reference scans exist, if so, get CC coefficients from there
    if MRSCont.flags.hasRef
        try
            if isSpecial
                % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
                cweights          = op_getcoilcombos(op_combinesubspecs(MRSCont.raw_ref_uncomb{kk}, 'diff'), 1, 'h');
            else
                cweights          = op_getcoilcombos(MRSCont.raw_ref_uncomb{kk},1,'h');
            end
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{ll,kk},1,'h',cweights);
            raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{ref_ll,kk},1,'h',cweights);
            if MRSCont.flags.isUnEdited
                raw_comb.flags.isUnEdited = 1;
                raw_ref_comb.flags.isUnEdited = 1;
            elseif MRSCont.flags.isMEGA
                raw_comb.flags.isMEGA = 1;
                raw_ref_comb.flags.isMEGA = 1;
            elseif MRSCont.flags.isHERMES
                raw_comb.flags.isHERMES = 1;
                raw_ref_comb.flags.isHERMES = 1;
            elseif MRSCont.flags.isHERCULES
                raw_comb.flags.isHERCULES = 1;
                raw_ref_comb.flags.isHERCULES = 1;
            end
            MRSCont.raw{ll,kk}     = raw_comb;
            MRSCont.raw_ref{ref_ll,kk} = raw_ref_comb;
            if MRSCont.raw_ref{ref_ll,kk}.subspecs > 1 && ~isSpecial && (length(size(MRSCont.raw_ref{ref_ll,kk}.fids)) > 2)
                if MRSCont.flags.isMEGA
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},2);
                    MRSCont.raw_ref{ref_ll,kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},2);
                    raw_ref_C               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},3);
                    raw_ref_D               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},4);
                    MRSCont.raw_ref{ref_ll,kk} = op_concatAverages(op_concatAverages(raw_ref_A,raw_ref_B),op_concatAverages(raw_ref_C,raw_ref_D));
                end
            else
                % Maintain spatial sub-spectra for SPECIAL-localized data.
                if ~isSpecial
                    MRSCont.raw_ref{kk}.subspecs = 1;
                    MRSCont.raw_ref{kk}.dims.subSpecs=0;
                end
            end
            
        catch
            
            % if wrong number of channels etc, use the metabolite scan itself
            if isSpecial
                % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
                cweights          = op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(MRSCont.raw_uncomb{ll,kk}), 'diff'), 0, 0.01, 2);
            else
                cweights          = op_getcoilcombos(MRSCont.raw_uncomb{ll,kk}, 1, 'h');
            end
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{ll,kk},1,'h',cweights);
            cweights            = op_getcoilcombos(MRSCont.raw_ref_uncomb{ref_ll,kk},1,'h');
            raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{ref_ll,kk},1,'h',cweights);
            if MRSCont.flags.isUnEdited
                raw_comb.flags.isUnEdited = 1;
                raw_ref_comb.flags.isUnEdited = 1;
            elseif MRSCont.flags.isMEGA
                raw_comb.flags.isMEGA = 1;
                raw_ref_comb.flags.isMEGA = 1;
            elseif MRSCont.flags.isHERMES
                raw_comb.flags.isHERMES = 1;
                raw_ref_comb.flags.isHERMES = 1;
            elseif MRSCont.flags.isHERCULES
                raw_comb.flags.isHERCULES = 1;
                raw_ref_comb.flags.isHERCULES = 1;
            end
            MRSCont.raw{kk}     = raw_comb;
            MRSCont.raw_ref{ref_ll,kk} = raw_ref_comb;
            if MRSCont.raw_ref{ref_ll,kk}.subspecs > 1 && ~isSpecial && (length(size(MRSCont.raw_ref{ref_ll,kk}.fids)) > 2)
                if MRSCont.flags.isMEGA
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},2);
                    MRSCont.raw_ref{ref_ll,kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},2);
                    raw_ref_C               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},3);
                    raw_ref_D               = op_takesubspec(MRSCont.raw_ref{ref_ll,kk},4);
                    MRSCont.raw_ref{ref_ll,kk} = op_concatAverages(op_concatAverages(raw_ref_A,raw_ref_B),op_concatAverages(raw_ref_C,raw_ref_D));
                end
            else
                % Maintain spatial sub-spectra for SPECIAL-localized data.
                if ~isSpecial
                    MRSCont.raw_ref{kk}.subspecs = 1;
                    MRSCont.raw_ref{kk}.dims.subSpecs=0;
                end
            end
        end
        
    else
        % if not, use the metabolite scan itself
        if isSpecial
            % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
            cweights          = op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(MRSCont.raw_uncomb{ll,kk}), 'diff'), 0, 0.01, 2);
        else
            cweights          = op_getcoilcombos(MRSCont.raw_uncomb{ll,kk}, 1, 'h');
        end
        raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{ll,kk},1,'h',cweights);
        if MRSCont.flags.isUnEdited
            raw_comb.flags.isUnEdited = 1;
        elseif MRSCont.flags.isMEGA
            raw_comb.flags.isMEGA = 1;
        elseif MRSCont.flags.isHERMES
            raw_comb.flags.isHERMES = 1;
        elseif MRSCont.flags.isHERCULES
            raw_comb.flags.isHERCULES = 1;
        end
        MRSCont.raw{ll,kk}     = raw_comb;
    end
    
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        if isSpecial
            % Workflow adopted from https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
            cweights_w          = op_getcoilcombos(op_combinesubspecs(MRSCont.raw_w_uncomb{w_ll,kk}, 'diff'), 1, 'h');
        else
            cweights_w          = op_getcoilcombos(MRSCont.raw_w_uncomb{w_ll,kk}, 1, 'h');
        end
        raw_w_comb          = op_addrcvrs(MRSCont.raw_w_uncomb{w_ll,kk},1,'h',cweights_w);
        raw_w_comb.flags.isUnEdited = 1;
        MRSCont.raw_w{w_ll,kk}   = raw_w_comb;
    end
end
%% Clean up and save

% Set flags
MRSCont.flags.coilsCombined     = 1;


end
