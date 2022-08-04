%% Project: Shape-based attitude planning %%
% Date: 31/01/22

%% Cost function %%
% Function to estimate the time of flight

% Inputs: - matrix I, the inertia dyadic of the system
%         - string cost, the cost function to be minimized
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - scalar mu, the gravitational parameter of the central body 
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - string basis, the polynomial basis to be used
%         - string dynamics, the dynamics formulation to be used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(cost, I, initial, final, n, tau, x, B, basis, dynamics)
    % Optimization variables 
    tf = x(end);                                        % Maneuver time

    % Minimize the cost function 
    switch (cost)
        case 'Minimum energy'
            % Minimize the control input
            P = reshape(x(1:end-1), [length(n), max(n)+1]);                 % Control points
            P = boundary_conditions(tf, n, initial, final, P, B, basis);    % Boundary conditions control points
            C = evaluate_state(P, B, n);                                    % State evolution
            u = acceleration_control(I, C, tf, dynamics);                   % Control input    
            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);                        % Dimensional control input
         
            r = trapz(tau,a)/tf;                                            % Minimum energy applied to the system

        case 'Minimum time'
            r = tf; 

        case 'Minimum power'
            % Minimize the control input
            P = reshape(x(1:end-1), [length(n), max(n)+1]);                 % Control points
            P = boundary_conditions(tf, n, initial, final, P, B, basis);    % Boundary conditions control points
            C = evaluate_state(P, B, n);                                    % State evolution
            [u, dv] = acceleration_control(I, C, tf, dynamics);             % Control input    
         
            r = abs(trapz(tau,dot(dv,u,1))/tf^2);                           % Minimum energy applied to the system

        otherwise
            error('No valid cost function was selected');
    end
end