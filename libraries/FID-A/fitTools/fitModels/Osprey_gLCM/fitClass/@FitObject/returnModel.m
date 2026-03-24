function model = returnModel(obj,step,secDim)
%%  returnModel(obj,step,secDim)
%   This method generates returns data and model spectra from the fit object.
%
%   USAGE:
%       obj.plotFit(step,secDim)
%
%   INPUTS:
%       step            = step to plot      
%       secDim          = spectrum along indirect dimension to plot % OPTIONS:   - [] default (plot all) 
%                                                                                - n (plot spectrum with index n)
%
%   OUTPUTS:
%       model struct
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%%  Diverge to default options if required 

    if nargin < 3
        secDim = 1;                                             % Set second dimensions to plot
        if nargin < 2
            step = obj.step;                                    % Set to last step
        end
    end

%%  Get fit data from object 

    data = fftshift(fft(obj.Data.fids,[],1),1);                                         % Get data matrix
    fit = obj.Model{step}.fit.fit;                                                      % Get fit matrix
    residual = obj.Model{step}.fit.residual;                                            % Get residual matrix           
    baseline = obj.Model{step}.fit.baseline;                                            % Get baseline matrix
    metabs = obj.Model{step}.fit.metabs;                                                % Get metabolite matrix

%% Get second dimension and fill struct
    if ~isnan(secDim)
        model.data = data(:,secDim);
        model.fit  = fit(:,secDim);
        model.residual = residual(:,secDim);
        model.baseline = baseline(:,secDim);
        model.metabs = metabs(:,:,secDim);
    else
        model.data = data(:,:);
        model.fit  = fit(:,:);
        model.residual = residual(:,:);
        model.baseline = baseline(:,:);
        model.metabs = metabs(:,:,:);
    end
end