%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, m0, Isp, T, tau, initial, final, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);
    tf = x(end-1);
    N = floor(x(end));

    % Equalities 
    ceq = [];

    % Boundary conditions points
    P(:,[1 2 end-1 end]) = boundary_conditions(mu, tf, n, initial, final, N, 'Orthogonal Bernstein');

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    u = acceleration_control(mu,C,tf);

    % Inequality (control authority)
    c = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-tf^2*T*ones(1,size(u,2));
end