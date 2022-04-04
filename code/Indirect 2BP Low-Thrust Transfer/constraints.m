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

    C = evaluate_state(P,B,n);

    % Boundary conditions
    ceq = [];
    ceq = [C(1:length(initial),1)-initial.';  ...
           C(1:length(final),end)-final.'];

    % Collocation equations do not apply to boundary conditions
%     C = C(:,2:end-1);
%     tau = tau(2:end-1);

    % Inequality (control authority)
    u = zeros(6,size(C,2));
    c = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-ones(1,size(u,2));

    % Relevant variables
    r = sqrt(C(1,:).^2+C(3,:).^2);
    m = m0-tf*Isp*tau;

    % Collocation equations
    f = [C(4,:); ...
         C(5,:)./C(1,:); ...
         C(6,:); ... 
         C(5,:).^2./C(1,:)-mu.*C(1,:)/r.^3; ... 
         -C(4,:).*C(5,:)./C(1,:); ...
         -mu.*C(3,:)./r.^3];

    A = zeros(size(f)).';

    D = C(13:18,:)-tf*(f+u./m);      % State dynamics
    L = C(19:24,:)+tf*A.';           % Co-state dynamics

    % Transversality constraints 
    T = [dot(C(7:12,1),f(:,1)); dot(C(7:12,end),f(:,end))];

    % Final constraints
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1]); reshape(L, [size(L,1)*size(L,2) 1]); reshape(T, [size(T,1)*size(T,2) 1])];
end