%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Initial fitting %%
% Function to estimate the trajectory approximation

% Inputs: - scalar n, the degree of the approximation 
%         - vector tau, the control parameter vector 
%         - C, the initial trajectory estimation
%         - string basis, the Bernstein polynomial basis to be used

% Outputs: - array P, the estimation of the boundary control points as a
%            cell
%          - array C, the initial estimation of the spacecraft state vector

function [P, C] = initial_fitting(n, tau, C, basis)
    % Preallocation of the control points and the polynomials
    P = zeros(length(n), max(n)+1); 
    B = state_basis(n, tau, basis);

    % Compute the position control points leveraging the complete state vector
    for i = 1:length(n)
        A = B{i}(1:n(i)+1,:);
        P(i,1:n(i)+1) = C(i,:)*pinv(A);
    end
    
    % Evaluate the state vector
    C = evaluate_state(P, B, n);
end