%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string dynamics, the independent variable parametrization to be
%           used
%         - string cost, the cost function to be minimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % Final time of flight 
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);
    C(6:10,:) = C(6:10,:)/tf;

    % Control input 
    [u, ~] = acceleration_control(mu, C, tf);
    u = u/tf^2;

    % Equalities 
    F = zeros(4,size(C,2));
    for i = 1:size(C,2)
        F(:,i) = KS_matrix(C(1:4,i))*C(6:9,i);
    end
    ceq = [mu*C(10,:)/4+dot(F(1:3,:),u(1:3,:),1)./F(4,:)];

    % Inequalities
    c = [dot(u(1:3,:),u(1:3,:),1)-(repmat(T,1,size(u,2))).^2];
end