%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string dynamics, the independent variable parametrization to be
%           used
%         - string cost, the cost function to be minimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-3), [length(n), max(n)+1]);     % Control points
    tf = x(end-2);                                      % Final time of flight 
    N = floor(x(end-1));                                % Optimal number of revolutions
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    [u, ~] = acceleration_control(mu, C, tf);

    % Equalities 
    ceq = [];

    c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*repmat(T,1,size(u,2)))];
end