function [sopt,vopt,s0]= kScaleOptimization_info_alpha(x,alpha,s0)
% Automatic tuning of the scale parameter for the exponentiated quadratic
% kernel: y = exp(d.^2/(2*s^2))
% FORMAT [sigma, value] = kScaleOptimization_info_alpha(d,alpha,s0)
% x     - distance among data points, usually: x = pdist2(x_in,x_in);
% alpha - Renyi's entropy alpha value
% s0    - starting point for the search
% sopt - achieved optimum value
% vopt - objective function value at sigma
%
% The tuning method used here is based on the maximization of information
% Potential variability as a function of the scale parameter, since
% lim_{s->0}{var{y(s)}} = 0 and lim_{s->inf}{var{y(s)}} = 0 and a
%__________________________________________________________________________
% Copyright (C) 2014 Signal Processing and Recognition Group
% David Cardenas & Andres Marino Alvarez Meza
% $Id: ScaleOptimization_info_alpha.m 2014-02-22 22:40:00 $

if nargin < 2
    s0 = median(x(:));
    alpha = 2;
elseif nargin < 3
    s0 = median(x(:));
end
    
f = @(s)obj_fun(s,x,alpha);
[sopt, vopt] = fminsearch(f,s0,alpha);
%[sopt, vopt] = fmincon(f,s0,[],[],[],[],eps,1e15);
if sopt<0
    sopt = abs(sopt);
end
 %%%%% objective function %%%%%%%%%%%
 function [v] = obj_fun(s,x,alpha)
 
 k = size(x,1)^(2-alpha)*(exp(-x.^2/(2*s^2))).^(alpha-1);
 v = - var(mean(k));