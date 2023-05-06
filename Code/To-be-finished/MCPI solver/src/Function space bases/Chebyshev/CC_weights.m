%% Project:  
% Sergio Cuevas del Valle
% Date: 06/11/22
% File: CC_weights
% Issue: 0 
% Validated: 

%% Clenshaw-Curtis weights %%
% This scripts provides the function to compute the Clenshaw-Curtis weights
% for a given domain interval and polynomial degree. Ref Greg von Winckel -
% 02/12/2005.

% Inputs: - vector x, the Legendre nodes at which the weights shall be
%           evaluated
%         - vector dP, the n-th Legendre derivative at the Legendre nodes

% Output: - vector w, the Legendre weights of interest

function [w, tau] = CC_weights(N)  
   % Constants 
   n = N-1; 

   % FFT preallocation of the CC weights and nodes
   c = zeros(N,2);
   c(1:2:N,1) = (2./[1 1-(2:2:n).^2]).';
   c(2,2) = 1; 
   f = real(ifft([c(1:N,:); c(n:-1:2,:)]));
   w = [f(1,1); 2*f(2:n,1); f(N,1)];
   tau = n*f(1:N,2);

   [tau,i] = sort(tau); 
   w = w(i);
end