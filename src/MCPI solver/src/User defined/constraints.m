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

function [c, ceq] = constraints(cost, mu, initial, final, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-8), [length(n), max(n)+1]);     % Control points for the control law
    lambda0 = x(end-7:end-2);                           % Initial co-state guesses
    tf = x(end-1);                                      % Final time of flight 
    T = x(end);                                         % Needed thrust vector

    % Trajectory evolution
    u = evaluate_state(P,B,n) / tf^2;
    
    % Integrate the trajectory
    tol = 1e-5;
    x0 = repmat(initial, size(u,2), 1);
    dyn = @(tau,x)dynamics(mu, tf, u, tau, x);
    [x, state] = MCPI(tau, x0, dyn, length(tau)-1, tol);

    % Equalities 
    % ceq = [x(end,[1 3])-final([1 3]) cos(x(end,2))-cos(final(2)) sin(x(end,2))-sin(final(2)) x(end,4:6)-final(4:6)];
    ceq = x(end,:)-final;

%     switch (cost)
%         case 'Minimum time'
%             % Integrate the co-states 
%             lambda = repmat(lambda0.', size(C,2), 1);
%             dynamics = @(tau, x)dual_dynamics(tf, C, x);
%             [lambda, ~] = MCPI(tau, lambda, dynamics, length(tau)-1, 1e-5);
%     
%             % Compute the Hamiltonian constraint
%             H = lambda(end,:)*C(4:9,end)+1;
% 
%         case 'Minimum fuel'    
%             % Integrate the co-states 
%             lambda = repmat(lambda0.', size(C,2), 1);
%             dynamics = @(tau, x)dual_dynamics(tf, C, x);
%             [lambda, ~] = MCPI(tau, lambda, dynamics, length(tau)-1, 1e-5);
% 
%             % Compute the primer vector constraint
%             Lambda = lambda(:,4:6)./sqrt(dot(lambda(:,4:6),lambda(:,4:6),2));
%             U = u./sqrt(dot(u,u,1));
%             H = 1+dot(Lambda, U.', 1);
% 
%         otherwise
%     end
% 
%     ceq = [ceq H];

    % Inequalities
    c = [u(1,:).^2+u(2,:).^2+u(3,:).^2-(repmat(T,1,size(u,2))).^2];
end

%% Auxiliary functions 
% Dual dynamics 
function [dlambda] = dual_dynamics(tf, C, lambda)
    % Preallocation 
    dlambda = zeros(size(C,2),6); 

    A = [zeros(3) eye(3); zeros(3,6)];

    for i = 1:size(C,2)
        r = norm(C(1:3,i));
        A(4:6,1:3) = [-1/r^3+3*C(1,i)^2/r^5 3*C(1,i)*C(2,i)/r^5 3*C(1,i)*C(3,i)/r^5; ...
                      3*C(1,i)*C(2,i)/r^5 -1/r^3+3*C(2,i)^2/r^5 3*C(2,i)*C(3,i)/r^5; ...
                      3*C(1,i)*C(3,i)/r^5 3*C(2,i)*C(3,i)/r^5 -1/r^3+3*C(3,i)^2/r^5]; 
        dlambda(i,:) = -tf*lambda(i,:)*A;
    end
end