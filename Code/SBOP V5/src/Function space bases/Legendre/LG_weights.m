%% Project:  
% Sergio Cuevas del Valle
% Date: 06/11/22
% File: LG_weights
% Issue: 0 
% Validated: 

%% Legendre-Gauss weights %%
% This scripts provides the function to compute the Legendre-Gauss weights
% for a given domain interval and polynomial degree 

% Inputs: - vector x, the Legendre nodes at which the weights shall be
%           evaluated
%         - vector dP, the n-th Legendre derivative at the Legendre nodes

% Output: - vector w, the Legendre weights of interest

function [w] = LG_weights(x, dP)    
    w = 2./((1-x.^2).*(dP.^2));
end