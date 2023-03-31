function [out, fs, phs, w, driftPre, driftPost] = op_probabSpecReg(in, seqType, echo, F0, noAlign)
% Spectral registration is a time-domain frequency-and-phase correction
% routine as per Near et al. (2015). Incorporates a multiplexed,
% probabilistic approach for aligning HERMES data (MM: 170609)
showPlots = 0;
if nargin <5
    noAlign = 0;
    if nargin <4
        F0 = nan;
        if nargin < 2
            echo = 1;    
        end
    end
end

% Check whether data is coil-combined. If not, throw error.
if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% Check whether data is averaged. If it is, throw warning and return unmodified data.
if in.dims.averages == 0 || in.averages < 2
    %DO NOTHING
    warning('No averages found.  Returning input without modification!');
    out = in;
    out.flags.freqcorrected = 1;
    fs  = 0;
    phs = 0;
    return
end

% Check whether input structure has sub-spectra.
% If it does, prepare to loop over all sub-spectra. Each sub-spectrum will
% get aligned separately.
if in.dims.subSpecs == 0
    numSubSpecs = 1;
    % Measure drift pre-alignment
    driftPre = op_measureDrift(in);
else
    numSubSpecs = in.sz(in.dims.subSpecs);
    % Measure drift pre-alignment
    for mm = 1:numSubSpecs
        in_sub = op_takesubspec(in, mm);
        driftPre{mm} = op_measureDrift(in_sub);
    end
end


% Pre-allocate memory.
fids    = zeros(in.sz(in.dims.t),1,numSubSpecs);
fs      = zeros(in.sz(in.dims.averages),numSubSpecs);
phs     = zeros(in.sz(in.dims.averages),numSubSpecs);
zMSE = zeros(numSubSpecs,in.sz(in.dims.averages));
CorrParsML = zeros(in.sz(in.dims.averages),2);
parsGuess = [0 0];

% Initialize common variables and loop over sub-spectra
t                 = in.t;
input.dwelltime   = in.dwelltime;

% Probability density function and parameter bounds
Cauchy = @(x,s,l) s./(pi.*(s.^2+(x-l).^2));
lb = [0 -Inf];
ub = [Inf Inf];

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Display','off','Algorithm','levenberg-marquardt');
% nlinopts = statset('nlinfit');
% nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-8,'TolFun',1e-8);
mleopts = statset('mlecustom');
mleopts = statset(mleopts,'MaxIter',400,'MaxFunEvals',800,'TolX',1e-6,'TolFun',1e-6,'TolBnd',1e-6);

% Set dimensions of figures of histograms
if showPlots == 1
    d.w = 0.6;
    d.h = 0.45;
    d.l = (1-d.w)/2;
    d.b = (1-d.h)/2;
end

for mm=1:numSubSpecs
    clear m;
    freq        = in.ppm;
    DataToAlign = in.fids(:,:,mm);
    % Use first n points of time-domain data, where n is the last point where SNR > 3
    signal = abs(DataToAlign(:,mm));
    noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
    SNR = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
    SNR = abs(diff(mean(SNR,2)));
    %SNR = abs(diff(mean(SNR,2))); This is how it's done in RobSpecReg
%     n = find(SNR > 1.5); % This is original Gannet
%     if isempty(n)
%         tMax = 100;
%     else
%         tMax = n(end);
%         if tMax < 100
%             tMax = 100;
%         end
%     end

    SNR = SNR(t <= 0.2);
    tMax = find(SNR > 1.5,1,'last');
    if isempty(tMax) || tMax < find(t <= 0.1,1,'last')
        tMax = find(t <= 0.1,1,'last');
    end

    % 'Flatten' complex data for use in nlinfit
    clear flatdata;
    flatdata(:,1,:) = real(DataToAlign(1:tMax,:));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,:));

    % Reference transient
    flattarget = median(flatdata,3); % median across transients
    target = flattarget(:);

    % Scalar to normalize transients (reduces optimization time)
    a = max(abs(target));
    % Pre-allocate memory
    if mm == 1
        params = zeros(size(flatdata,3), 2);
        MSE = zeros(1, size(flatdata,3));
    end

     % Frame-by-frame determination of frequency of residual water and Cr (if HERMES or GSH editing)
    if strcmp(seqType, 'HERMES') || strcmp(seqType, 'HERCULES')
        F0freqRange = freq - 3.02 >= -0.15 & freq - 3.02 <= 0.15;
    else
        F0freqRange = freq - 3.22 >= -0.15 & freq - 3.22 <= 0.15;
    end
    [~,FrameMaxPos] = max(abs(real(in.specs(F0freqRange,:))),[],1);
    F0freqRange = freq(F0freqRange);
    F0freq = F0freqRange(FrameMaxPos);

    % Starting values for optimization
    if isnan(F0)

        f0 = F0freq * in.txfrq * 1e-6;
        f0 = f0 - f0(1);
    else
        f0 = F0;
        f0 = f0 - f0(1);
    end
    phi0 = zeros(size(f0));
    x0 = [f0(:) phi0(:)];

    % Determine frequency and phase offsets by spectral registration
     if noAlign == 0
        time = 0:input.dwelltime:(length(target)/2 - 1)*input.dwelltime;
        reverseStr = '';
        for jj = 1:size(flatdata,3)
            if echo
                msg = sprintf('\nSpectral registration - Fitting transient: %d', corrloop);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end

            transient = squeeze(flatdata(:,:,jj));

            fun = @(x) SpecReg(transient(:)/a, target/a, time, x);
            params(jj,:) = lsqnonlin(fun, x0(jj,:), [], [], lsqnonlinopts);

            f   = params(jj,1);
            phi = params(jj,2);
            m_c = complex(flatdata(:,1,jj), flatdata(:,2,jj));
            m_c = m_c .* exp(1i*pi*(time'*f*2+phi/180));
            m(:,jj) = [real(m_c); imag(m_c)];
            resid = target - m(:,jj);
            MSE(jj) = sum(resid.^2) / (length(resid) - 2);

%             x0(jj,2) = phi; % Update phase value according to the last iteration

    %         [parsFit(corrloop,:), ~, ~, ~, MSE(corrloop)] = nlinfit(input, target, @FreqPhaseShiftNest, parsGuess, nlinopts);
    %         parsGuess = parsFit(corrloop,:);
        end
         % Prepare the frequency and phase corrections to be returned by this
        % function for further use
        fs(:,mm)  = params(:,1);
        phs(:,mm) = params(:,2);
    end

    % Probability distribution of frequency offsets (estimated by maximum likelihood)
    start = [iqr(fs(:,mm))/2, median(fs(:,mm))];
    [fs_p(:,mm), fs_p_ci(:,:,mm)] = ...
        mle(fs(:,mm), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    fs_fx(:,mm) = ...
        linspace(1.5*min(fs(:,mm)), 1.5*max(fs(:,mm)), 1e3);
    fs_pdf(:,mm) = Cauchy(fs_fx(:,mm), fs_p(1,mm), fs_p(2,mm));

    % Probability distribution of phase offsets (estimated by maximum likelihood)
    start = [iqr(phs(:,mm))/2, median(phs(:,mm))];
    [phs_p(:,mm), phs_p_ci(:,:,mm)] = ...
        mle(phs(:,mm), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    phs_fx(:,mm) = ...
        linspace(1.5*min(phs(:,mm)), 1.5*max(phs(:,mm)), 1e3);
    phs_pdf(:,mm) = Cauchy(phs_fx(:,mm), phs_p(1,mm), phs_p(2,mm));

             
    zMSE = zscore(MSE); % standardized MSEs
    

    % Apply frequency and phase corrections
     for jj = 1:size(flatdata,3)
        % Default correction
%         fids(:,jj,mm) = in.fids(:,jj,mm) .* ...
%             exp(1i*params(jj,1)*2*pi*t') * exp(1i*pi/180*params(jj,2));

        % Freq/phase correction + Cauchy pdf location parameter shift
        fids(:,jj,mm) = in.fids(:,jj,mm) .* ...
            exp(1i*(params(jj,1)-fs_p(2))*2*pi*t') * exp(1i*pi/180*(params(jj,2)-phs_p(2)));
     end
end

if echo
fprintf('\n');
end

%%% 2. CALCULATE THE DRIFT AFTER CORRECTIONS
out_temp=in;
out_temp.fids=fids;
%re-calculate Specs using fft
out_temp.specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);
% Measure drift post-alignment
if in.dims.subSpecs == 0
    numSubSpecs = 1;
    driftPost = op_measureDrift(out_temp);
else
    numSubSpecs = in.sz(in.dims.subSpecs);
    % Measure drift post-alignment
    for mm = 1:numSubSpecs
        out_temp_sub = op_takesubspec(out_temp, mm);
        driftPost{mm} = op_measureDrift(out_temp_sub);
    end
end


%%% 3. RE-RUN THE RANKING ALGORITHM TO DETERMINE WEIGHTS %%%

% Determine optimal iteration order by calculating a similarity metric (mean squared error)
w = cell(1,numSubSpecs);

for mm=1:numSubSpecs
    DataToWeight = fids(:,:,mm);
    D = zeros(size(DataToWeight,2));
    for jj = 1:size(DataToWeight,2)
        for kk = 1:size(DataToWeight,2)
            tmp = sum((real(DataToWeight(1:tMax,jj)) - real(DataToWeight(1:tMax,kk))).^2) / tMax;
            if tmp == 0
                D(jj,kk) = NaN;
            else
                D(jj,kk) = tmp;
            end
        end
    end
    d = nanmedian(D);
    if isnan(sum(d))

        d(isnan(d)) = 1;
    end
    w{mm} = 1./d.^2;
    w{mm} = w{mm}/sum(w{mm});
    w{mm} = repmat(w{mm}, [size(DataToWeight,1) 1]);
        
    
    % Apply the weighting
fids_out(:,:,mm) = w{mm} .* DataToWeight;

end

% Or remove based on zMSE score
% for mm=1:numSubSpecs
%     DataToWeight = fids(:,:,mm);
%     w{mm} = zMSE < 3;
%     w{mm} = w{mm}/sum(w{mm});
%     w{mm} = repmat(w{mm}, [size(DataToWeight,1) 1]);
%         
%     % Apply the weighting
%     fids_out(:,:,mm) = w{mm} .* DataToWeight;
% end

% Perform weighted averaging
% No need for 'mean', since the weights vector w is normalized
fids = squeeze(sum(fids_out, in.dims.averages));

%%% 3. WRITE THE NEW STRUCTURE %%%
%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

%change the dims variables.
if in.dims.t>in.dims.averages
    dims.t=in.dims.t-1;
else
    dims.t=in.dims.t;
end
if in.dims.coils>in.dims.averages
    dims.coils=in.dims.coils-1;
else
    dims.coils=in.dims.coils;
end
dims.averages=0;
if in.dims.subSpecs>in.dims.averages
    dims.subSpecs=in.dims.subSpecs-1;
else
    dims.subSpecs=in.dims.subSpecs;
end
if in.dims.extras>in.dims.averages
    dims.extras=in.dims.extras-1;
else
    dims.extras=in.dims.extras;
end

%re-calculate the sz variable
sz=size(fids);

% Add NIfTI-MRS provenance
% Generate fields for provenance
fields.Method   = 'Frequency and phase correction';
fields.Details  = ['Probabilistic multi-step registration to median (Mikkelsen et al. 2018), dim = DIM_DYN'];
in = op_add_analysis_provenance(in,fields);
fields.Method   = 'Signal averaging';
fields.Details  = ['Similarity metric weighted sum (Mikkelsen et al. 2020), dim = DIM_DYN'];
in = op_add_analysis_provenance(in,fields);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;
out.averages=1;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqcorrected=1;
out.flags.averaged=1;



end

function out = SpecReg(data, target, time, x)

f   = x(1);
phi = x(2);

data = reshape(data, [length(data)/2 2]);
fid = complex(data(:,1), data(:,2));

fid = fid .* exp(1i*pi*(time'*f*2+phi/180));
fid = [real(fid); imag(fid)];

out = target - fid;

end

