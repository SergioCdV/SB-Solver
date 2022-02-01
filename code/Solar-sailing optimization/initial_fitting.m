%% Project: 
% Date: 01/02/22

%% Initial fitting %%
% Function to estimate the trajectory approximation

% Inputs: - scalar n, the degree of the approximation 
%         - vector tau, the control parameter vector 
%         - Capp, the initial trajectory estimation

% Outputs: - array B, the Bernstein polynomials basis in use
%          - array P, the estimation of the boundary control points
%          - array C, the initial estimation of the spacecraft state vector

function [B, P, C] = initial_fitting(n, tau, C)
    % The Bernstein-basis polinomials for the increasd order are calculated
    B = [bernstein_basis(n,tau); bernstein_derivative(n,tau,1); bernstein_derivative(n,tau,2)];

    % State vector fitting
    P = C(1:3,:)*pinv(B(1:n+1,:));
    C = repmat(P,1,3)*B;
end