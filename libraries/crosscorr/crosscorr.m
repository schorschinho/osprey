function xcorr = crosscorr(x,y)
%% xcorr = crosscorr(x,y)
%   This function calculates the cross-correlation between two vectors x and y.
%
%   USAGE:
%       xcorr = crosscorr(x,y)
%
%   INPUTS:
%       x     = a vector
%       y     = another vector
%
%   OUTPUTS:
%       xcorr = cross correlation between vector x and y.
%
%   AUTHOR:
%       Dr. Helge Zöllner (Johns Hopkins University, 2020-03-15)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on the follwoing code:
%       Sidhanta Kumar Panda (2020). 
%       Cross-correlation (https://www.mathworks.com/matlabcentral/fileexchange/43681-cross-correlation)
%       MATLAB Central File Exchange. Retrieved May 15, 2020.
%
%   HISTORY:
%       2020-05-15: First version of the code.
%% Calculate cross-correlation
%%% 1.PREPARE DATA %%%
h=fliplr(y);
lx=length(x);
lh=length(h);
n=lx+lh-1;
hh=[h zeros(1,n-lh)];
xx=zeros(n);
xx(1:lx,1)=x;

%%% 2.CALCULATE CROSS-CORRELATION %%%
for i=2:n
    for j=2:n
        xx(j,i)=xx(j-1,i-1);
      
    end
end
xcorr=xx*hh';

end
