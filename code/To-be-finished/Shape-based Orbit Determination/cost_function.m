%% Project: 
% Date: 31/01/22

%% Cost function %%
% Function to estimate the time of flight

% Inputs: - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - scalar mu, the gravitational parameter of the central body 
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - string basis, the polynomial basis to be used
%         - string method, the parameter distribution to be used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(mu, n, t, r, tau, x, B, basis, method)
    % Minimize the control input
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % The Bézier control points

    % Boundary conditions
    P(:,[1 end]) = boundary_conditions(n, t, r, basis);

    % State evolution
    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(mu,C,tf,method); 

    % Control cost
    a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);            % Dimensional acceleration
    r = trapz(tau,a);                                   % Cost function
end