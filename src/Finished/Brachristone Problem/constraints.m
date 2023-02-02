%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(g, tf, tau, initial, final, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-m-1), [3, max(n)+1]);
    C = evaluate_state(P,B,n);
    u = reshape(x(end-m:end-1), [1 m]);

    % Non-linear inequality
    c = [];

    % Boundary conditions
    ceq = [];
    C(4:6,:) = C(4:6,:);
    ceq = [C(1:3,1)-initial(1:3).'; C(1:2,end)-final(1:2).'];

    % Dynamic constraints  
    C = C(:,2:end-1);
    u = u(:,2:end-1);
    D = [C(4:5,:)-x(end)*C(3,:).*[sin(u); -cos(u)]; C(6,:)-x(end)*g*cos(u)];

    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end