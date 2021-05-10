function map = osp_generateMRSImap(MRSCont, kk, nominator_spec,denominator_spec, nominator, denominator)
%% map = osp_generateMRSImap(MRSCont, kk, nominator_spec,denominator_spec, nominator, denominator, figTitle)
%   Creates a matrix with the corresponding metabolite map
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%                                 'mm' re_mm
%
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided

if nargin<6
    denominator = [];
    if nargin<5
        nominator = {'Cr', 'PCr'}; 
        if nargin<4
            denominator_spec = 'none'; 
            if nargin < 3
                nominator_spec = 'off';
                if nargin < 2
                    kk = 1;
                    if nargin<1
                        error('ERROR: no input Osprey container specified.  Aborting!!');
                    end
                end
            end
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
map_dims = size(MRSCont.quantify.amplMets{kk}.(nominator_spec).(nominator{1}));
nominator_map = zeros([1,map_dims]);
denominator_map = zeros([1,map_dims]);
if  (MRSCont.flags.isMRSI == 1)
    for nom = 1 : length(nominator)
        nominator_map(1,:,:) = nominator_map(1,:,:) + reshape(MRSCont.quantify.amplMets{kk}.(nominator_spec).(nominator{nom}),[1 map_dims]);
    end
    if ~isempty(denominator)
        for denom = 1 : length(denominator)
            denominator_map(1,:,:) = denominator_map(1,:,:) + reshape(MRSCont.quantify.amplMets{kk}.(denominator_spec).(denominator{denom}),[1 map_dims]);
        end
        map = nominator_map ./ denominator_map;
        map(denominator_map == 0) = 0;
    else
        map = nominator_map;
    end    
end

map = permute(map, [2 3 1]);
map = rot90(map,2);
map = flipud(map);
end

   