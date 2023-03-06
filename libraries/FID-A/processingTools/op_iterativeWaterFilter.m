%op_iterativeWaterFilter.m
%Georg Oeltzschner, Johns Hopkins University 2020.
%
% USAGE:
% out=op_iterativeWaterFilter(out,wlim,Kinit,M,plot_bool);
% 
% DESCRIPTION:
% This function is a wrapper to attempt error-free water removal using
% multiple instances of the op_removeWater function. This is necessary
% because at any given number of frequency components, the HSVD filter
% might end up removing the entire FID instead.
%
% To circumvent this problem, we begin with an initial guess of 'Kinit'
% components. Should the FID end up empty after filtering, this function
% starts at 1.5*Kinit components, and decreases the Kinit variable by one
% until it arrives at a non-empty FID.

% At each step, the water signal is being removed using HSVD method 
% described by H. BARKHUIJSEN et al. 1987.
% 
% INPUTS:
% in        = MRS data structure used by FID-a toolkit. Data should be
%             pre-processed, for example by out = run_pressproc(filename)
% wlim      = This is the frequency limits of the water peak to be fitted in
%             ppm. (default = [4.4 5]
% Kinit     = The number of frequency components in the data model This parameter
%             might have to be played with. (default is 20).
% M         = M is the integer number of columns in the henkel matrix. Note: L
%             is the number of rows and L+M=N where N is the number of data
%             points. For best results 0.5<=L/M<=2. (default M= .75*length.
% plot_bool = if 1, water fit is plotted (default =1)
%
% OUTPUTS:
% out       = New spectrum without the water peak in the as a FID-A structure

function [out] = op_iterativeWaterFilter(in,wlim,Kinit,M,plot_bool)

% set default values ( intended for seimens data with 4096 data points)
if nargin<5
    plot_bool =0;
    if nargin<4
        M = floor(in.sz(1)*.75);
        if nargin<3
            Kinit=30;   % this value might have to be played with
            if nargin<2
                wlim= [4.4 5];
            end
        end
    end
end

out_temp = cell(in.subspecs,1);
for ss = 1 : in.subspecs
    if in.subspecs == 1
        out_temp{ss} = in;
    else
        out_temp{ss} = op_takesubspec(in,ss);
    end
    % Make the first attempt with the initial guess for the number of
    % components to be removed
    [out_temp{ss},~,~]   = op_removeWater(out_temp{ss},wlim,Kinit,M,plot_bool); % Remove the residual water

    % If the resulting FID is empty...
    if isnan(real(out_temp{ss}.fids))
        % ... increase the number of components to 1.5 times the initial guess
        K_rr = 1.5*Kinit;
        % Run again
        while (isnan(real(out_temp{ss}.fids(1))) && (K_rr > 0))
            if in.subspecs == 1
                out_temp{ss} = in;
            else
                out_temp{ss} = op_takesubspec(in,ss);
            end
            [out_temp{ss},~,~]   = op_removeWater(out_temp{ss},wlim,K_rr,M,plot_bool); % Remove the residual water
            % If the resulting FID is still empty, decrease the number of
            % components and try again until it is not empty.
            K_rr = K_rr-1;
        end
    end

    if isnan(real(out_temp{ss}.fids))
        out_temp{ss}.fids = out_temp{ss}.fids;
        out_temp{ss}.specs = out_temp{ss}.specs;
    end
    out_temp{ss}     = op_fddccorr(out_temp{ss},100); % Correct back to baseline   
end
% Use the first non-empty FID to return
out     = out_temp{1};
watersupp = out_temp{1}.watersupp;
out = rmfield(out,'watersupp');
out.watersupp{1} = watersupp;
for ss = 2 : in.subspecs
    out = op_mergesubspec(out,out_temp{ss});
    out.watersupp{ss} = out_temp{ss}.watersupp;
end
if exist('K_rr')
    K = K_rr;
else
    K = Kinit;
end
fields.Method   = 'Nuisance peak removal';
fields.Details  = ['HSVD removal of residual water (Barkhuijsen et al 1987) using ' num2str(K) ' components'];
out = op_add_analysis_provenance(out,fields);



end

