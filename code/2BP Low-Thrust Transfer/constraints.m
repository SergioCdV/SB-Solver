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
    P = reshape(x(1:end-1-2*m), [4, max(n)+1]);
    u = reshape(x(end-2*m:end-1), [2 m]);
    tf = x(end);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    ceq = [C([1 3 4],1)-initial([1 3 4]).';  ...
           C([1 3 4],end)-final([1 3 4]).'];

    % Inequality (control authority)
    c = sqrt(u(1,:).^2+u(2,:).^2)-T*ones(1,size(u,2));
    c = [];

    % Dynamic constraints
    m = m0-tf*Isp*tau;
    D = [C(5,:)+tf*( -C(3,:)); ...
         C(7,:)+tf*( -C(4,:).^2./C(1,:)+mu./C(1,:).^2-u(1,:)./m ); ...
         C(8,:)+tf*(  C(3,:).*C(4,:)./C(1,:)-u(2,:)./m ); ... 
         C(6,:)+tf*( -C(4,:)./C(1,:) )]; 

    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end