%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 18/08/22
% File: mfls_optimization.m 
% Issue: 0 
% Validated: 18/08/22

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar selector, to constrain the needed acceleration vector
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(wp, wv, kap, phi, psi, tau, T, initial, final, n, x, B, basis)
    % Extract the optimization variables
    P = reshape(x(1:end-1), [length(n), max(n)+1]);     % Control points
    tf = x(end);                                        % Final time of flight 

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    [u] = acceleration_control(wp, wv, kap, phi, psi, C, tf, tau);

    % Inequality (control authority)
    a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Dimensional acceleration

    ceq = [];
    c = a-T*ones(1,size(u,2));
end