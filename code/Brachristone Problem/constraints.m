%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(g, tf, tau, initial, final, n, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-1), [2, max(n)+1]);
    C = evaluate_state(P,B,n);

    % Non-linear inequality
    c = [];

    % Boundary conditions
    ceq = [C(1:2,1)-initial(1:2).'; norm(C(3:4,1)); C(1:2,end)-final(1:2).'];

    C = C(:,2:end-1);

    % Dimensionalising 
    C(3:4,:) = C(3:4,:) / (tf * x(end));
    C(5:6,:) = C(5:6,:) / (tf * x(end))^2;

    % Dynamic constraints  
    u = atan2(C(3,:),-C(4,:));
    v = sqrt(C(3,:).^2+C(4,:).^2);
    du = -(C(5,:).*C(4,:)-C(6,:).*C(3,:))./v.^2; 

    D = C(5:6,:)-v.*du.*[cos(u);sin(u)]-g*cos(u).*[sin(u);-cos(u)];
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end