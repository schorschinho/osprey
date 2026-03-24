% osp_spiral.m
% Helge Zollner, Johns Hopkins University 2020.
%
% USAGE:
% [x, y] = osp_spiral(n);
% 
% DESCRIPTION:
% This function calcualtes the x,y coordinates of the nth point of a ulam
% spiral. This is usful if you want to fit MRSI data following a spiral
% trajectory.
%
% 
% INPUTS:
% n           = nth point of the spiral

% OUTPUTS:
% x           = Cartesian x coordinate
% y           = Cartesian y coordinate
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2021-01-07)
%       hzoelln2@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-01-07: First version of the code.

function[x,y] = osp_spiral(n)
        k=ceil((sqrt(n)-1)/2);
        t=2*k+1;
        m=t^2 ;
        t=t-1;
        if n>=m-t 
            x= k-(m-n);
            y=-k;        
        else
            m=m-t ;            
            if n>=m-t 
                x=-k;
                y=-k+(m-n);
            else
                m=m-t; 
                if n>=m-t 
                    x=-k+(m-n);
                    y=k; 
                else
                    x=k;
                    y=k-(m-n-t); 
                end
            end
        end        
end