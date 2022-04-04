%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [r] = cost_function(T, x, B, m, n, time)
    % Maximize the orbital radius transfer
    P = reshape(x(1:end-1), [length(n), max(n)+1]);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    % Minimize the Hamiltonian free-time constraint
    r = 0;
end