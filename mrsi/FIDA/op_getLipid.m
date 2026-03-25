%op_getLipid.m
%Jay Hennessy, McGill University 2017.
%
% USAGE:
% out=op_getLipid(out,wlim,Kinit,M,plot_bool);
% 
% DESCRIPTION:
% This function genrates lipid basis from MRS data using HSVD method 
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
% K         = The number of frequency components used to fit the data.

function [fid_lip, spec_lip ] = op_getLipid(in,liplim,Kinit,M)

% set default values ( intended for seimens data with 4096 data points)
if nargin<4
    M = floor(in.sz(1)*.75);
    %M = 1500;  % good for data sets of 2048 points
    %M = 3000;  % good for data sets of 4096 points
    if nargin<3
        Kinit=30;   % this value might have to be played with
        if nargin<2
            liplim= [0.3 1.9];
        end
    end
end

% set parameters
N =in.sz(1);
dt = in.dwelltime;
t=in.t(1:N);
ppm = in.ppm;

fid = in.fids(1:N)';
spec = in.specs;

% Make a hankel data matrix
H = hankel(fid(1:M),fid(M:end));

% calculate the svd
[U,S,V] = svd(H);

% initialize while loop
amp = 0;
count = 0;

    
%initialize K
K = Kinit;

% truncate the data
Uk = U(:,1:K);
Sk = S(1:K,1:K);
Vk = V(:,1:K);

% get he eigenvalues of the transform matrix
Utk = Uk(2:end,:);
Ubk = Uk(1:end-1,:);
Eh = Utk\Ubk;

Eeigs = eig(Eh');

% convert eigenvalues to poles to get freq and damping factor
[lip_model,alpha_model] = cart2pol(real(Eeigs), imag(Eeigs));
lip = lip_model/dt;
freqs = lip/(2*pi);
alpha = (alpha_model-1)/dt;

% make a model guess using only damping factor and freq
fid_temp = exp((-alpha + (1i*lip))*t);  % this is mostly zero except for 1 very small row

% do a least square fit of your model guess to the data
phamp = fid_temp'\fid';
% NofNan = size(find(phamp==0))


% convert the eigenvalues to phase and amplitude of the model
[ph,amp] = cart2pol(real(phamp), imag(phamp));


% use alpitude, phase, damping factor and frequency to model data
fid_components = (amp.*exp(1i*ph))'*exp((-alpha + (1i*lip))*t);
spec_components=fftshift(fft(fid_components',[],1),1);



% Model the lipid signal 
lipppm = -(freqs)/(in.txfrq/1000000)+in.centerFreq;
lipid = find(lipppm>liplim(1) & lipppm<liplim(2)); % find the frequencies components associated with water

fid_lip = (amp(lipid).*exp(1i*ph(lipid)))'*exp((-alpha(lipid) + (1i*lip(lipid)))*t);


%create water signal to be stored
spec_lip=fftshift(fft(fid_lip',[],1),1);

end

