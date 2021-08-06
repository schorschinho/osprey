function [CRLBstruct] = fit_OspreyNoLS_CRLB(basisSet,fitParams);


nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions

ampl = fitParams.ampl;

gaussLBCol = 2 * ones(1,nMets + nMM);
freqShiftCol = 3 * ones(1,nMets + nMM);
asymCol = freqShiftCol(end) + 1;
lorentzLBCol = asymCol(end) + 1 : (asymCol(1) + nMets + nMM);
freqShiftColmets = lorentzLBCol(end) + 1 : (lorentzLBCol(end) + nMets + nMM);
AmplCol = freqShiftColmets(end) + 1 : (freqShiftColmets(end) + nMets + nMM);
SplineAmplCol = AmplCol(end) + 1 : (AmplCol(end) + length(fitParams.beta_j));

D = fitParams.J;
sigma = fitParams.stdNoise;
D = D * sigma *sigma;


%calculate the fisher matrix
fisher = (1./(sigma^2)) * real(D.'*D);

CRLB   = sqrt(pinv(fisher));
CRLB = diag(CRLB);

%Calculate for common metabolite combinations
D_cut = D;
comb_signal = 0;
cut = [];

% Create combos
% tNAA  NAA + NAAG
idx_1 = find(strcmp(basisSet.name,'NAA'));
idx_2 = find(strcmp(basisSet.name,'NAAG'));
if  ~isempty(idx_1) && ~isempty(idx_2)
    idx_3 = find(strcmp(basisSet.name,'tNAA'));
    if isempty(idx_3)
        basisSet.name{length(basisSet.name)+1} = 'tNAA';
    end
    D_comb = D(:,AmplCol(idx_1)) +D(:,AmplCol(idx_2));
    cut(end+1) = AmplCol(idx_1);
    cut(end+1) = AmplCol(idx_2);
    D_cut(:,end+1) = D_comb;
    ampl(end+1) = ampl(idx_1) + ampl(idx_2);
    comb_signal = comb_signal + 1;
end
% Glx Glu+Gln
idx_1 = find(strcmp(basisSet.name,'Glu'));
idx_2 = find(strcmp(basisSet.name,'Gln'));   
if  ~isempty(idx_1) && ~isempty(idx_2)
    idx_3 = find(strcmp(basisSet.name,'Glx'));
    if isempty(idx_3)
        basisSet.name{length(basisSet.name)+1} = 'Glx';
    end
    D_comb = D(:,AmplCol(idx_1)) +D(:,AmplCol(idx_2));
    cut(end+1) = AmplCol(idx_1);
    cut(end+1) = AmplCol(idx_2);
    D_cut(:,end+1) = D_comb;
    ampl(end+1) = ampl(idx_1) + ampl(idx_2);
    comb_signal = comb_signal + 1;
end
% tCho GPC+PCh
idx_1 = find(strcmp(basisSet.name,'GPC'));
idx_2 = find(strcmp(basisSet.name,'PCh'));
if  ~isempty(idx_1) && ~isempty(idx_2)
    idx_3 = find(strcmp(basisSet.name,'tCho'));
    if isempty(idx_3)
        basisSet.name{length(basisSet.name)+1} = 'tCho';
    end
    D_comb = D(:,AmplCol(idx_1)) +D(:,AmplCol(idx_2));
    cut(end+1) = AmplCol(idx_1);
    cut(end+1) = AmplCol(idx_2);
    D_cut(:,end+1) = D_comb;
    ampl(end+1) = ampl(idx_1) + ampl(idx_2);
    comb_signal = comb_signal + 1;
end
% tCr Cr+PCr
idx_1 = find(strcmp(basisSet.name,'Cr'));
idx_2 = find(strcmp(basisSet.name,'PCr'));
if  ~isempty(idx_1) && ~isempty(idx_2)
    idx_3 = find(strcmp(basisSet.name,'tCr'));
    if isempty(idx_3)
        basisSet.name{length(basisSet.name)+1} = 'tCr';
    end
    D_comb = D(:,AmplCol(idx_1)) +D(:,AmplCol(idx_2));
    cut(end+1) = AmplCol(idx_1);
    cut(end+1) = AmplCol(idx_2);
    D_cut(:,end+1) = D_comb;
    ampl(end+1) = ampl(idx_1) + ampl(idx_2);
    comb_signal = comb_signal + 1;
end

%PE+EA
idx_1 = find(strcmp(basisSet.name,'PE'));
idx_2 = find(strcmp(basisSet.name,'EA'));
if  ~isempty(idx_1) && ~isempty(idx_2)
    idx_3 = find(strcmp(basisSet.name,'tEA'));
    if isempty(idx_3)
        basisSet.name{length(basisSet.name)+1} = 'tEA';
    end
    D_comb = D(:,AmplCol(idx_1)) +D(:,AmplCol(idx_2));
    cut(end+1) = AmplCol(idx_1);
    cut(end+1) = AmplCol(idx_2);
    D_cut(:,end+1) = D_comb;
    ampl(end+1) = ampl(idx_1) + ampl(idx_2);
    comb_signal = comb_signal + 1;
end

% %GABA+coMM3
% if strcmp(MRSCont.opts.fit.coMM3, '1to1GABA') % fixed GABA coMM3 model      
%     idx_1 = find(strcmp(basisSet.name,'GABA'));
%     if  ~isempty(idx_1)
%         idx_3 = find(strcmp(basisSet.name,'GABAplus'));
%         if isempty(idx_3)
%             basisSet.name{length(basisSet.name)+1} = 'GABAplus';
%         end
%         idx_GABAp = find(strcmp(basisSet.name,'GABAplus'));
%         GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1,:);
%         MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp,:) = GABAp;
%     end    
% else if strcmp(MRSCont.opts.fit.coMM3, '3to2MM') % fixed MM09 coMM3 model
%         idx_1 = find(strcmp(basisSet.name,'GABA'));
%         idx_2 = find(strcmp(basisSet.name,'MM09'));
%         if  ~isempty(idx_1) && ~isempty(idx_2)
%             idx_3 = find(strcmp(basisSet.name,'GABAplus'));
%             if isempty(idx_3)
%                 basisSet.name{length(basisSet.name)+1} = 'GABAplus';
%             end
%             idx_GABAp = find(strcmp(basisSet.name,'GABAplus'));
%             GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1,:) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2,:);
%             MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp,:) = GABAp;
%         end          
% else % Models with a separate comMM3 function or without a co-edited MM function         
%         idx_1 = find(strcmp(basisSet.name,'GABA'));
%         idx_2 = find(strcmp(basisSet.name,'MM3co'));
%         if  ~isempty(idx_1) && ~isempty(idx_2)
%             idx_3 = find(strcmp(basisSet.name,'GABAplus'));
%             if isempty(idx_3)
%                 basisSet.name{length(basisSet.name)+1} = 'GABAplus';
%             end
%             idx_GABAp = find(strcmp(basisSet.name,'GABAplus'));
%             GABAp = MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_1,:) + MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_2,:);
%             MRSCont.quantify.amplMets{kk}.(getResults{ll})(idx_GABAp,:) = GABAp;
%         end
% end
% end

D_cut(:,cut) = [];

fisher_cut = (1./(sigma^2)) * real(D_cut.'*D_cut);
CRLB_cut   = sqrt(pinv(fisher_cut));
CRLB_cut = diag(CRLB_cut);
n_CRLB = length(CRLB_cut);

if comb_signal > 0
    CRLB(end+1 : end + comb_signal) = CRLB_cut(n_CRLB - comb_signal + 1 : n_CRLB);
end

CRLB = real(CRLB);

% Create a struct with all CRLB values 
% Absolute values
CRLBstruct.abs.ph0 = CRLB(1);
CRLBstruct.abs.GaussLB = CRLB(gaussLBCol(1));
CRLBstruct.abs.freqShift = CRLB(freqShiftCol(1));
CRLBstruct.abs.asym = CRLB(asymCol);
CRLBstruct.abs.lorentzLB = CRLB(lorentzLBCol);
CRLBstruct.abs.freqShiftmets = CRLB(freqShiftColmets);
CRLBstruct.abs.amps = [CRLB(AmplCol); CRLB(end - comb_signal + 1 : end)];
CRLBstruct.abs.splineAmps = CRLB(SplineAmplCol);

% Relative values
CRLBstruct.rel.ph0 = 100./fitParams.ph0 .* CRLBstruct.abs.ph0;
CRLBstruct.rel.GaussLB = 100/ fitParams.gaussLB .* CRLBstruct.abs.GaussLB;
CRLBstruct.rel.freqShift = 100./ fitParams.freqShift .* CRLBstruct.abs.freqShift;
CRLBstruct.rel.asym = 100./ fitParams.asym .* CRLBstruct.abs.asym;
CRLBstruct.rel.lorentzLB = 100./  fitParams.lorentzLB .* CRLBstruct.abs.lorentzLB;
CRLBstruct.rel.freqShiftmets = 100./ fitParams.freqShiftmets .* CRLBstruct.abs.freqShiftmets;
CRLBstruct.rel.amps = 100./ ampl .* CRLBstruct.abs.amps;
CRLBstruct.rel.splineAmps = 100./ fitParams.beta_j .* CRLBstruct.abs.splineAmps;

end