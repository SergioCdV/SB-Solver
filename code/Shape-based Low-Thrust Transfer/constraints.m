%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string method, the parameter distribution to be used

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, T, initial, final, n, x, B, basis, method, tau)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % Final time of flight 
    N = floor(x(end));                                  % Optimal number of revolutions

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Radius constraints
    r = sqrt(C(1,:).^2+C(3,:).^2);
    r0 = sqrt(initial(1)^2+initial(3)^2);
    rf = sqrt(final(1)^2+final(3)^2);

    % Control input 
    u = acceleration_control(mu,C,tf,method);

    % Equalities 
    ceq = [];

    % Inequality (control authority)
    c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-tf^2*T*ones(1,size(u,2)) r-2*max([r0 rf])];
end