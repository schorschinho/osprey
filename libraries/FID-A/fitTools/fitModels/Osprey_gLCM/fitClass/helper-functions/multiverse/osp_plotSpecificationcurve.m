function[out] = osp_plotSpecificationcurve(MRSCont, Metab, referencing, Dataset, SubSpec, Ex)
%% function[out] = osp_plotSpecificationcurve(MRSCont, Metab, referencing, Dataset, SubSpec, Ex)
%
% Description: Generates a specification curve analysis using an Osprey
%              container. The models are ordered according to
%              "referencing", and arranged according to model procedure
%              filenames (must be BIDS-formatted, i.e. key1-val1_key2-val2)
%
% Input:     MRSCont = Osprey container
%            Metab = specific metabolite on which to perform SCA (default = 'tNAA')
%            referencing = Referencing mode (default = 'tCr')
%            Dataset = The index of the spectrum (or spectra) to view
%            SubSpec = Sub spectrum
%            Ex = Extra dim
%
% Example usage:
%                   osp_plotSpecificationcurve(MRSCont, 'GABA', 'tCr', 2);
%                   (plots the GABA/tCr ratio from the GABA-edited difference spectrum)
%                   
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
MRSCont = [];
Metab  {mustBeText} = 'tNAA';
referencing {mustBeText} = 'amplMets';
Dataset = 0;
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

AllBasisFuncs = unique([MRSCont.quantify.names.metab{1:end}]);
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
                    if isinf(Val)
                        Amp(jj_BasisFunc, jj_part, jj_models) = 1000;
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

if Dataset==0
    % We take the median across the spectra to provide a single metric for SCA
    Amp_Partmed = squeeze(median(Amp,2,'omitnan'));
else
    % We take the median across the specified datasets (just the value, if one dataset specified)
    Amp_Partmed = squeeze(median(Amp(:,Dataset,:),2,'omitnan'));
end

MetabBool = matches(AllBasisFuncs,Metab);
if ~any(MetabBool)
    error('No metabolite amplitudes found labeled with "%s"',Metab)
end
[Amp_Partmed, SrtInd] = sort(Amp_Partmed(MetabBool,:));

%% Create an informative struct assuming BIDS-formatted key-value pairs in the model JSON filenames
% Filename format:  key1-value-1_key2-value2_..._keyN-value-N.json

[~,List] = fileparts(MRSCont.opts.fit.ModelProcedure.metab(:,SubSpec));
for jj_keys=1:length(List)
    Pairs = strsplit(List{jj_keys},'_');
    for KK=1:length(Pairs)
        Singles = strsplit(Pairs{KK},'-');
        if ~(length(Singles)==2)
            error('Filename %i incompatible formating:\n\t%s',jj_keys,List{jj_keys});
        else
            KVP(jj_keys).(Singles{1}) = Singles{2};
        end
    end
end

%% Plot the specification curve

if MRSCont.flags.isGUI
    out = figure( 'Visible', 'off' );
else
    out = figure;
end

% First plot the ordered amplitdues across all models
subplot(3,1,1)
plot(Amp_Partmed,'ko-','LineWidth',1.5,'MarkerFaceColor',[0 0 0],'MarkerSize',12)
set(gca,'XTick',[])
XL = 1-0.2-length(SrtInd)/10;
xlim([XL length(SrtInd)+0.2])
ylabel([Metab, ' [',referencing ,']'])
set(gca,'FontSize',11);


% Then visualize the model key-value pairs
subplot(3,1,2:3)
Model_Keys = fieldnames(KVP);
Cols = cbrewer('qual','Dark2',max([length(Model_Keys),3]),'linear');

KeySep = 0;
ValSep = 0;
for jj_keys=1:length(Model_Keys)
    Model_Vals = unique({KVP.(Model_Keys{jj_keys})});
    for jj_vals = 1:length(Model_Vals)
        KeyVal_bool = matches({KVP.(Model_Keys{jj_keys})}, Model_Vals{jj_vals});

        XVal = find(KeyVal_bool(SrtInd));
        YVal = KeySep + ValSep;

        plot(XVal, YVal,'o','color',Cols(KeySep+1,:),'MarkerFaceColor',Cols(KeySep+1,:),'MarkerSize',20);
        hold on
        text(XL,YVal,['  ' Model_Vals{jj_vals}],'FontSize',11,'HorizontalAlignment','left','color',Cols(KeySep+1,:)) % Model Value
        plot([1,NModels], [YVal,YVal],'k') % Add horizontal guide line
        
        ValSep = ValSep+1; % Gap between sets of model vals
    end
    text(XL,YVal-(length(Model_Vals)-1)/2,Model_Keys{jj_keys},'FontSize',25,'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90,'color',Cols(KeySep+1,:)) % Model key

    KeySep = KeySep+2; % Gap between sets of model keys
end

xlim([XL NModels+0.2])
ylim([-1 YVal+1])
set(gca, "YTick", [], "XTick", [])

end