function [MRSCont] = osp_combineCoils(MRSCont,kk)
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


%% Calculate coil combination weights
if nargin<2
% Loop over all datasets
    for kk = 1:MRSCont.nDatasets
        % Check if reference scans exist, if so, get CC coefficients from there
        if MRSCont.flags.hasRef
            cweights            = op_getcoilcombos(MRSCont.raw_ref_uncomb{kk},1,'h');
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{kk},1,'h',cweights);
            raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{kk},1,'h',cweights);
            MRSCont.raw{kk}     = raw_comb;
            MRSCont.raw_ref{kk} = raw_ref_comb;
            if MRSCont.raw_ref{kk}.subspecs > 1
                if MRSCont.flags.isMEGA
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    MRSCont.raw_ref{kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    raw_ref_C               = op_takesubspec(MRSCont.raw_ref{kk},3);
                    raw_ref_D               = op_takesubspec(MRSCont.raw_ref{kk},4);                    
                    MRSCont.raw_ref{kk} = op_concatAverages(raw_ref_A,raw_ref_B,raw_ref_C,raw_ref_D);
                end
            end            
        else
            % if not, use the metabolite scan itself
            cweights            = op_getcoilcombos(MRSCont.raw_uncomb{kk},1,'h');
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{kk},1,'h',cweights);
            MRSCont.raw{kk}     = raw_comb;
        end

        % Now do the same for the (short-TE) water signal
        if MRSCont.flags.hasWater
            cweights_w          = op_getcoilcombos(MRSCont.raw_w_uncomb{kk},1,'h');
            raw_w_comb          = op_addrcvrs(MRSCont.raw_w_uncomb{kk},1,'h',cweights_w);
            MRSCont.raw_w{kk}   = raw_w_comb;
        end
    end
else
    % Check if reference scans exist, if so, get CC coefficients from there
    if MRSCont.flags.hasRef
        try
            cweights            = op_getcoilcombos(MRSCont.raw_ref_uncomb{kk},1,'h');
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{kk},1,'h',cweights);
            raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{kk},1,'h',cweights);
            MRSCont.raw{kk}     = raw_comb;
            MRSCont.raw_ref{kk} = raw_ref_comb;
            if MRSCont.raw_ref{kk}.subspecs > 1
                if MRSCont.flags.isMEGA
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    MRSCont.raw_ref{kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    raw_ref_C               = op_takesubspec(MRSCont.raw_ref{kk},3);
                    raw_ref_D               = op_takesubspec(MRSCont.raw_ref{kk},4);                    
                    MRSCont.raw_ref{kk} = op_concatAverages(op_concatAverages(raw_ref_A,raw_ref_B),op_concatAverages(raw_ref_C,raw_ref_D));
                end
            end
        catch
        % if wrong number of channels etc, use the metabolite scan itself
            cweights            = op_getcoilcombos(MRSCont.raw_uncomb{kk},1,'h');
            raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{kk},1,'h',cweights);
            MRSCont.raw{kk}     = raw_comb;
            cweights            = op_getcoilcombos(MRSCont.raw_ref_uncomb{kk},1,'h');
            raw_ref_comb        = op_addrcvrs(MRSCont.raw_ref_uncomb{kk},1,'h',cweights);
            MRSCont.raw_ref{kk} = raw_ref_comb;    
            if MRSCont.raw_ref{kk}.subspecs > 1
                if MRSCont.flags.isMEGA
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    MRSCont.raw_ref{kk} = op_concatAverages(raw_ref_A,raw_ref_B);
                else
                    raw_ref_A               = op_takesubspec(MRSCont.raw_ref{kk},1);
                    raw_ref_B               = op_takesubspec(MRSCont.raw_ref{kk},2);
                    raw_ref_C               = op_takesubspec(MRSCont.raw_ref{kk},3);
                    raw_ref_D               = op_takesubspec(MRSCont.raw_ref{kk},4);                    
                    MRSCont.raw_ref{kk} = op_concatAverages(op_concatAverages(raw_ref_A,raw_ref_B),op_concatAverages(raw_ref_C,raw_ref_D));
                end
            end
        end
    else
    % if not, use the metabolite scan itself
        cweights            = op_getcoilcombos(MRSCont.raw_uncomb{kk},1,'h');
        raw_comb            = op_addrcvrs(MRSCont.raw_uncomb{kk},1,'h',cweights);
        MRSCont.raw{kk}     = raw_comb;
    end

    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        cweights_w          = op_getcoilcombos(MRSCont.raw_w_uncomb{kk},1,'h');
        raw_w_comb          = op_addrcvrs(MRSCont.raw_w_uncomb{kk},1,'h',cweights_w);
        MRSCont.raw_w{kk}   = raw_w_comb;
    end
end    
%% Clean up and save

% Set flags
MRSCont.flags.coilsCombined     = 1;
% Delete un-combined data to free up memory
raw_fields = {'raw_uncomb','raw_ref_uncomb','raw_w_uncomb'};
for kk = 1:length(raw_fields)
    if isfield(MRSCont, raw_fields{kk})
        MRSCont = rmfield(MRSCont, raw_fields{kk});
    end
end

    
end