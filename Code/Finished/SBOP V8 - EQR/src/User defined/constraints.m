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

function [c, ceq] = constraints(mu, initial, final, tf, time_free, uprev, B, basis, n, L, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    thetaf = x(end-1);                                  % Final insertion phase
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    P = boundary_conditions(n, initial(1:5), final(1:5), P, B, basis);

    % Compute the longitude evolution 
    theta0 = initial(end)-thetaf*L(1);
    L = theta0+thetaf*L;

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    [u, ~] = acceleration_control(mu, C, L, uprev);
 
    % Kinematic constraint 
    w = 1+C(2,:).*cos(L)+C(3,:).*sin(L);
    res = sqrt(mu*C(1,:)).*(w./C(1,:)).^2;
    for i = 1:size(C,2)
        B = control_input(mu, C(:,i)); 
        res(i) = 1/(res(i)+B(6,3)*u(3,i));
    end

    % Inequalities
    c = [u(1,:).^2+u(2,:).^2+u(3,:).^2-(thetaf*repmat(T,1,size(u,2))).^2 -res];

    % Equalities
    ceq = [cos(L(end))-cos(final(6)) sin(L(end))-sin(final(6))];

    if (~time_free)
        ceq = [ceq tf-thetaf*trapz(tau,res)];
    end
end