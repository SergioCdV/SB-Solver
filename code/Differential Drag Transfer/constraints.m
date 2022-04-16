%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
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

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, n, x, B, basis)
    % Extract the optimization variables
    P = reshape(x(1:end-1), [length(n), max(n)+1]);     % Control points
    tf = x(end);                                        % Final time of flight 

    % Equalities 
    ceq = [];

    % Boundary conditions points
    P(:,[1 2 end-1 end]) = boundary_conditions(tf, n, initial, final, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    u = acceleration_control(mu,C,tf);

    % Inequality (control authority)
    c = [];
end