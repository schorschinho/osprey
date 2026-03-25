function out = op_echo_filter(in, R2_L, R2p_L, sigma)
% APPLYLINEBROADENING Apply mixed Lorentzian and Gaussian line broadening to spin echo FID
%
% Inputs:
%   FID    - Original FID signal (complex array)
%
%   R2_L   - Lorentzian R2 rate (irreversible, in s^-1)
%   R2p_L  - Lorentzian R2' rate (reversible, in s^-1)
%   sigma  - Gaussian broadening parameter (in s^-1)
%
% Output:
%   FID_broadened - FID with applied line broadening
%
% Example:
%   tau = 0.070;  % 70 ms (TE = 140 ms)
%   dt = 0.001;   % 1 ms dwell time
%   t = -0.065:dt:0.260;  % Start 65ms before echo, end 260ms after
%   FID = exp(1i*2*pi*100*t);  % 100 Hz oscillation
%   FID_broad = applyLinebroadening(FID, t, tau, 3, 2, 2.5);

out = in;

% Convert relative time to absolute time
% t is relative to TE, so absolute time = t + TE
tau = in.te/2/1000;

t_abs = in.t + in.te/1000;

% Calculate phase evolution for spin echo
% Phase accumulates until tau, then reverses until TE, then accumulates again
phi = zeros(size(in.t));

for i = 1:length(in.t)
    if t_abs(i) < tau
        % Before 180° pulse: phase accumulates
        phi(i) = t_abs(i);
    elseif t_abs(i) < in.te/1000
        % Between 180° pulse and echo: phase reverses
        phi(i) = 2*tau - t_abs(i);
    else
        % After echo: phase accumulates again
        phi(i) = t_abs(i) - in.te/1000;
    end
end

% Calculate decay factors
% R2_L decay (irreversible - accumulates throughout)
decay_R2_L = exp(-R2_L * t_abs)';

% R2'_L decay (reversible - follows phase evolution)
decay_R2p_L = exp(-R2p_L * abs(phi))';

% Gaussian decay (reversible - follows phase evolution)
decay_Gaussian = exp(-0.5 * (sigma * phi).^2)';

% Apply all decay factors to the FID
fids = in.fids .* decay_R2_L .* decay_R2p_L .* decay_Gaussian;

specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

out.fids=fids;
out.specs=specs;


end