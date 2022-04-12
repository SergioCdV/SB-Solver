%% Project: 
% Date: 31/01/22

%% Cost function %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - scalar m, the number of grid points
%         - vector n, containing the order of approximation of each phase
%           space
%         - vector tau, the collocation grid points

% Outputs: - scalar r, the cost function to be optimized

function [r] = cost_function(x, B, m, n, tau)
    % Control points and control input
    u = reshape(x(end-3*m:end-1), [3 m]);
    tf = x(end);

    % Cost functino
    r = trapz(tau, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)/tf);
end