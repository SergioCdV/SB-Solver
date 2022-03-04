%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, m0, Isp, T, tf, tau, initial, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-m), [2, max(n)+1]);
    theta = reshape(x(end-m+1:end), [m, 1]).';

    % Non-linear inequality
    c = [];

    % Boundary conditions
    C = evaluate_state(P,B,n);
    C(3:4,:) = C(3:4,:)/tf;
    C(5:6,:) = C(3:4,:)/tf^2;
    ceq = [C(1:4,1)-initial(1:4).'; C(3:4,end)-[0; sqrt(mu/C(1,end))]];

    % Dynamic constraints  
    D = [C(5,:)-C(4,:).^2./C(1,:)+mu./C(1,:).^2-(T*sin(theta))./(m0-Isp*tf*tau); C(6,:)+C(3,:).*C(4,:)./C(1,:)-(T*cos(theta))./(m0-Isp*tf*tau)];
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end