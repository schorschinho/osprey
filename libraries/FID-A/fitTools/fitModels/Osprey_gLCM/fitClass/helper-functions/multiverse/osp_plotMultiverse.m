function[out] = osp_plotMultiverse(MRSCont, referencing, PartInd, SubSpec, Ex)
%% function[out] = osp_plotMultiverse(MRSCont, referencing, PartInd, SubSpec, Ex)
%
% Description: Plot distributions of model parameters (for e.g., amplitudes
% and CRLBS) across several model definitions. 
%
% Input:     MRSCont = MRS container with multiple rows of model
%               procedures, which has been run up to OspreyQuantify.
% Optional:  referencing = PArameter to use (e.g. amplMets for metabolites,
%               tCr for creatine ratios, or CRLB.
%            PartInd = The index for a particular dataset, or a vector, to
%               average across several, or 0 to average across all.
%            SubSpec = The sub spectrum of interest
%            Ex = The extra dimension of interest
% Output:    out = figure handle
%
% Example usage:
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
MRSCont = [];
referencing {mustBeText} = 'amplMets';
PartInd = 0;
SubSpec {mustBeInteger} = 1;
Ex {mustBeInteger} = 1;
end

%% Handle errors and set up
if size(MRSCont.opts.fit.ModelProcedure.metab,1)==1
    error('You need to specify more than 1 model JSON in dim 1!')
elseif ~(MRSCont.flags.didFit)    
    error('Must run OspreyFit first!')
elseif ~(MRSCont.flags.didQuantify)
    error('Must run OspreyQuantify first!')
end

NModels = size(MRSCont.opts.fit.ModelProcedure.metab,1);
NPart = MRSCont.nDatasets(1);

%% Retrieve the amplitudes across models into a 3D array
AllBasisFuncs = unique([MRSCont.quantify.names.metab{1,SubSpec,:,:}]);
if matches(referencing,'tCr')
    AllBasisFuncs(matches(AllBasisFuncs,{'tCr','Cr','PCr'})) = []; % Remove Cr measures if we're investigating tCr ratios
end

Amp = zeros(length(AllBasisFuncs),NPart,NModels).*nan;
for jj_BasisFunc=1:length(AllBasisFuncs)
    for jj_part=1:NPart
        for jj_models=1:NModels
            Bool = matches(MRSCont.quantify.names.metab{1,SubSpec,Ex,jj_models}, AllBasisFuncs{jj_BasisFunc});
            if any(Bool)
                if matches(referencing,'amplMets')
                    Amp(jj_BasisFunc, jj_part, jj_models) = MRSCont.quantify.amplMets{1,jj_part,SubSpec,Ex,jj_models}.metab(Bool);
                elseif matches(referencing, 'CRLBs')
                    Val = MRSCont.quantify.CRLB{1,jj_part,SubSpec,Ex,jj_models}.metab(Bool);
                    if isinf(Val) % If CRLB is infinite (zero amplitude), then replace it with nan for plotting
                        Amp(jj_BasisFunc, jj_part, jj_models) = nan;
                    else
                        Amp(jj_BasisFunc, jj_part, jj_models) = Val;
                    end
                else
                    Amp(jj_BasisFunc, jj_part, jj_models) = MRSCont.quantify.metab.(referencing){1,jj_part,SubSpec,Ex,jj_models}(Bool);
                end
            end 
        end
    end
end

if PartInd == 0
    Amp = squeeze(median(Amp,2));
else
    Amp = squeeze(median(Amp(:,PartInd,:),2));
end

[~,Ind] = sort(median(Amp,2,'omitnan'),'descend'); % plot in decending order

Amp_srt = Amp(Ind,:).';
Labels = AllBasisFuncs(Ind);

if MRSCont.flags.isGUI
    out = figure( 'Visible', 'off' );
else
    out = figure;
end

JitteredBoxPlot(Amp_srt,[],'Xlabels',Labels) % Boxplot of measures

% Labels
if matches(referencing, 'CRLBs')
    YL = 'CRLBs [%]';
    ylim([0 100]);
    yline(20)
else
    YL = ['Amp. [',referencing ,']'];
end
ylabel(YL);

set(gca,'FontSize',11);

end


function[BC, Output] = JitteredBoxPlot(InMat, Colors, varargin)
%% function[BC] = JitteredBoxPlot(InMat, Colors)
%
% Description: Plot boxplots alongside scattered individual points. (Taken
% from CWDAVIESJENKINS/CWDJ_MiscUtilities
%
% Input:     InMat = Matrix with columns representing seperate measures or
%               a cell array of vectors
% Optional:  Colors = RGB color vectors, with rows corresponding to the
%               columns of InMat (default uses cbrewer)
%            varargin = Additional (pairwise options)
% Output:    BC = Boxchar array
%            Output = A struct containing the median, quartiles, outliers,
%               and whisker locations for each input dataset.


%% Check inputs

% Convert to cell array (if matrix or table)
if iscell(InMat)
    S = length(InMat);
else
    S = size(InMat);
    M=InMat;InMat = cell(1,S(2));
    for JJ=1:S(2)
        Vec= M(:,JJ);
        Vec = Vec(~isnan(Vec));
        InMat{JJ} = Vec;
    end
end

% Create/check color matrix
if ~exist('Colors','var') || isempty(Colors)
    Colors = cbrewer('qual','Dark2',length(InMat),'linear'); % Use cbrewer function to generate RGB vectors using the "Dark2" set
else
    ColSize = size(Colors);
    if ColSize(1)<length(InMat)
        error('Color vector is too small (%i, compared to %i matrix columns)',ColSize,S(1));
    end
end

%% Manage Varargin

%%%%% Inititalize default parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Connection = [];                        % Array defining group connections
ConnectionAlpha = 0.4;                  % Alpha of the connecting lines
Width = 0.2;                            % Width of the scatter
Offset = 0.35;                          % How shifted the scatter is (set to 0 to plot over)
PlotOutliers = false;                   % Bool—whether to to plot outliers in boxchart
PointAlpha = 0.5;                       % Alpha of the scatter points
PointSize = 25;                         % Size of the scatter points
PointShape = repmat({'o'},S(2),1);      % Shape of scatter point (cell array)
Xlabels = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overwrite default parameters, if specified in varargin
if nargin>2
    L = length(varargin)/2;
    Args = reshape(varargin,2,L)';
    for JJ=1:L
        switch Args{JJ,1}
            case "Connection"
                Connection = Args{JJ,2};
            case "ConnectionAlpha"
                ConnectionAlpha = Args{JJ,2};
            case "Width"
                Width = Args{JJ,2};
            case "Offset"
                Offset = Args{JJ,2};
            case "PlotOutliers"
                PlotOutliers = Args{JJ,2};
            case "PointAlpha"
                PointAlpha = Args{JJ,2};
            case "PointSize"
                PointSize = Args{JJ,2};
            case "PointShape"
                PointShape = Args{JJ,2};
            case "Xlabels"
                Xlabels = Args{JJ,2};
            otherwise
                warning('Unknown option: %s',Args{JJ,1})
        end
    end
end

%% Loop to plot box charts and scatter

for JJ=1:length(InMat)
    % Populate secondary output struct with relevant descriptive statistics
    Output(JJ).Median = median(InMat{JJ});
    Output(JJ).Quartile = [prctile(InMat{JJ},25), prctile(InMat{JJ},25)];
    Output(JJ).IsOutlier = isoutlier(InMat{JJ},'quartiles');
    Output(JJ).Whisker = [min(InMat{JJ}(~Output(JJ).IsOutlier)), max(InMat{JJ}(~Output(JJ).IsOutlier))];

    if PlotOutliers
        bc(JJ) = boxchart(JJ*ones(1,length(InMat{JJ})),InMat{JJ},'BoxFaceColor',Colors(JJ,:),'MarkerColor',Colors(JJ,:),'BoxWidth',Width);
    else
        bc(JJ) = boxchart(JJ*ones(1,length(InMat{JJ})),InMat{JJ},'BoxFaceColor',Colors(JJ,:),'BoxFaceColor',Colors(JJ,:),'BoxWidth',Width,'MarkerStyle','none');
    end
    hold on
    
    if PointSize>0

        Jitter{JJ} = rand([1,length(InMat{JJ})]).*Width - Width/2 + JJ-Offset;

        scatter(Jitter{JJ}, InMat{JJ},'MarkerEdgeColor',Colors(JJ,:),...
                                      'MarkerFaceColor',Colors(JJ,:),...
                                      'MarkerEdgeAlpha',PointAlpha,...
                                      'MarkerFaceAlpha',PointAlpha,...
                                      'Marker',PointShape{JJ},...
                                      'SizeData',PointSize);
    end
end
Ax = gca;
Ax.XAxis.Visible = 'off'; % remove x-axis


%% Add optional xlabels

if ~isempty(Xlabels)
    for JJ=1:length(InMat)
        text(JJ,0,[Xlabels{JJ},'  '],'rotation',65,'horizontalalignment','right','fontsize',12,'color',Colors(JJ,:));
    end
end

%% Plot connected points

if ~isempty(Connection)
    Sc = size(Connection);
    for JJ=1:Sc(1)
        L1 = length(Jitter{Connection(JJ,1)});
        L2 = length(Jitter{Connection(JJ,2)});
        if ~(L1 == L2)
            error("Can't plot this connection! Entry %i has different lengths: %i and %i",JJ,L1,L2);
        end
        for KK=1:L1
            plot([Jitter{Connection(JJ,1)}(KK), Jitter{Connection(JJ,2)}(KK)],[InMat{Connection(JJ,1)}(KK), InMat{Connection(JJ,2)}(KK)],'Color',[0,0,0,ConnectionAlpha],'LineWidth',0.2)
        end
    end
end

% Return figure handle only if outputs are requested
if nargout > 0
    BC = bc;
end

end