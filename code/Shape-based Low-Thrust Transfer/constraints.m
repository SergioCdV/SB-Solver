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
    P = reshape(x(1:end-1), [length(n), max(n)+1]);
    tf = x(end);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    r = sqrt(C(1,:).^2+C(3,:).^2);
    m = m0-tf*Isp*tau;

    % Control input 
    u = [m.*( -C(10,:)+tf*( C(5,:).^2./C(1,:)-mu.*C(1,:)/r.^3 ) ); ...
         m.*( -C(11,:)+tf*( -C(4,:).*C(5,:)./C(1,:) ) ); ... 
         m.*( -C(12,:)+tf*( -mu.*C(3,:)/r.^3 ) )];

    ceq = [];
    ceq = [C(1:length(initial),1)-initial.';  ...
           C(1:length(final),end)-final.'];

    % Inequality (control authority)
    c = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-ones(1,size(u,2));

    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end