%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - mu, the gravitational parameter of the system 
%         - scalar T, the maximum allowed acceleration
%         - vector initial, the initial 6 by 1 state vector (heliocentric position and velocity)
%         - vector final, the final 6 by 1 state vector (heliocentric position and velocity)
%         - vector n, the approximation degree to each position coordinate
%         - scalar m, the number of grid points of interest
%         - vector x, containing the optimization degrees of freedom
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, T, initial, final, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-1-3*m), [length(n), max(n)+1]);
    u = reshape(x(end-3*m:end-1), [3 m]);
    tf = x(end);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    ceq = [];
    ceq = [C(1:length(initial),1)-initial.';  ...
           C(1:length(final),end)-final.'];

    % Inequality (control authority)
    c = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-T*ones(1,size(u,2));

    % Dynamic constraints
    r = sqrt(C(1,:).^2+C(3,:).^2);
    D = [C(7,:)+tf*( -C(4,:) ); ...
         C(8,:)+tf*( -C(5,:)/C(1,:) ); ...
         C(9,:)+tf*( -C(6,:) ); ...
         C(10,:)+tf*( -C(5,:).^2./C(1,:)+mu.*C(1,:)/r.^3-u(1,:) ); ...
         C(11,:)+tf*(  C(4,:).*C(5,:)./C(1,:)-u(2,:) ); ... 
         C(12,:)+tf*(  mu.*C(3,:)/r.^3-u(3,:) )];

    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end