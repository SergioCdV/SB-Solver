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
%         - matrix D, the differentiation matrix
%         - vector m, the number of nodes
%         - vector tau, the normalized independent variable
%         - vector x, the degree of freedom to be optimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, D, m, tau, x)
    % Extract the optimization variables
    C = reshape(x(1:6*(m+1)), [6 length(tau)]);                   % State points
    u = reshape(x(end-3*(m+1)-1:end-2), [3 length(tau)]);         % Control points
    tf = x(end-1);                                                % Final time of flight
    T = x(end);                                                   % Needed thrust vector

    % Boundary equalities 
    bcs = [C(1:6,1)-initial; C([1 3:6],end)-final([1 3:6]); cos(C(2,end))-cos(final(2)); sin(C(2,end))-sin(final(2)) ];

    % Path constraints 
    h = u(1,:).^2+u(2,:).^2+u(3,:).^2-(repmat(T,1,size(u,2))).^2;

    % Dynamics 
    res = acceleration_control(mu, D, C, u, tf, tau); 
    res = reshape(res, [], 1);

    % Final constraints
    ceq = [bcs; res];              % Boundary and dynamics constraints
    c = h;                      % Path inequalities
end