%op_SeparateMaxEcho.m
%
% USAGE:
% out=op_SeparateMaxEcho(in,side,tstart);
% 
% DESCRIPTION:
% 
% INPUTS:
% in        = MRS data structure used by FID-a toolkit. Data should be
%             pre-processed, for example by out = run_pressproc(filename)
% side      = which side of the echo do you want. See vatiations below
% tstart    = When did the ADC open (in ms).
%
%
% OUTPUTS:
% out       = New spectrum from left or right side of echo

function [out ] = op_SeparateMaxEcho(in,side,tstart)

 % Calculate time to before echo
tleft = in.te-tstart;
tleft = tleft/1000;
pointsToPick = tleft/in.dwelltime;

out = in;

switch side
    case 'left'   
        fids = zeros(in.sz);
        % fids(1:pointsToPick,:,:,:) = conj(flip(in.fids(1:pointsToPick,:,:,:),1));
        fids(1:pointsToPick,:,:,:) = ((in.fids(1:pointsToPick,:,:,:)));
    case 'flipleft'   
        fids = zeros(in.sz);
        fids(1:pointsToPick,:,:,:) = conj(flip(in.fids(1:pointsToPick,:,:,:),1));
    case 'right'
        fids = zeros(in.sz);
        fids(1:in.sz(1)-pointsToPick+1,:,:,:) = in.fids(pointsToPick:end,:,:,:);
    case 'add'
        fids = zeros(in.sz);
        fids(1:in.sz(1)-pointsToPick+1,:,:,:) = in.fids(pointsToPick:end,:,:,:);
        % fids(2:pointsToPick,:,:,:) = fids(2:pointsToPick,:,:,:) + conj(flip(in.fids(2:pointsToPick,:,:,:),1));
        fids(2:pointsToPick,:,:,:) = fids(2:pointsToPick,:,:,:) + ((in.fids(2:pointsToPick,:,:,:)));
    case 'shift'
        sz = in.sz;
        sz(1) = 2*in.sz(1)-pointsToPick;
        fids = zeros(sz);
        fids(in.sz(1)-pointsToPick+1:sz(1),:,:,:) = in.fids(:,:,:,:);
        out.sz = sz;
        f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):...
            in.spectralwidth/(sz(1)):...
            (in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];

        out.ppm=f/(out.Bo*42.577);
        out.ppm=out.ppm+out.centerFreq;
        out.t=[0:out.dwelltime:(sz(1)-1)*out.dwelltime];
end


out.fids = fids;
out.specs=fftshift(fft(out.fids,[],1),1);

end

