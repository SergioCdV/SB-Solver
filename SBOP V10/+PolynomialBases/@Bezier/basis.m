%% Project: 
% Date: 29/01/2022

%% Bernstein basis
% Function for calculation the Bernstein polynomial basis of degree n

% Inputs:  - scalar n, the cardinal of the basis
%          - vector tvec, the control parameter t vector 

% Output: - array B, the Bernstein polynomial basis of dimensions    
%           n x length(t)

% Iterative portion of De Casteljau's algorithm
function [B] = basis(order, u)    
    % Find number of steps (time increments)
    steps = length(u);
    
    % Initialize variable for n-order curve
    B = zeros(order+1,steps); 
    
    % Calculation of bezier curve, recursively
    for i = 0:order
        % Indexing variable 
        j = i+1; 
        
        % binomial coefficient (n i)
        bin = factorial(order)/(factorial(i)*factorial(order-i));
        
        % Bernstein basis polynomial
        B(j,:) = bin.*u.^i.*(1-u).^(order-i);
    end
end