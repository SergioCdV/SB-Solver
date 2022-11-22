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
    P = reshape(x(1:end-3), [length(n), max(n)+1]);     % Control points
    tf = x(end-2);                                      % Final time of flight 
    thetaf = x(end-1);                                  % Final insertion phase
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, thetaf, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    [u, ~] = acceleration_control(mu, C, tf);

    % Equalities 
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*cos(C(6,:));
    f = tf*[zeros(2,size(C,2)); sqrt(mu*C(1,:)).*(w./C(1,:)).^2];
    res = zeros(3,size(C,2)); 
    for i = 1:size(C,2)
        B = control_input(mu, C(:,i)); 
        res(:,i) = C([7 11 12],i)/tf-B([1 5 6],:)*u(:,i)-f(:,i);
    end

    res = dot(res,res,1); 
    ceq = [cos(thetaf)-cos(final(6)) sin(thetaf)-sin(final(6)) res];

    % Inequalities
    c = [u(1,:).^2+u(2,:).^2+u(3,:).^2-(repmat(T,1,size(u,2))).^2];
end