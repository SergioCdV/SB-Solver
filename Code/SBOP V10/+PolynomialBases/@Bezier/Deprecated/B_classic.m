%% Project: 
% Date: 29/01/2022

%% B_classic
% This function uses the classical implementation for calculating a Bezier curve as defined by the determined points.
% The time vector represents the number of points at which the curve B will be calculated

% Inputs: - array points, the Bézier curve control points 
%         - vector tvec, the control parameter t vector 

% Output: - the resulting Bézier curve B

function [B] = B_classic(points,tvec)
    % Determine order
    n = length(points)-1;

    % Generate the polynomial basis 
    bern = bernstein_basis(n,tvec);
    
    % Bézier curve as a dot product
    B = points*bern;
end