%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 
%         - string basis, the polynomial basis to be used
%         - string dynamics, the independent variable parametrization to be
%           used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(cost, mu, initial, final, n, tau, x, B, basis, dynamics)
    % Optimization variables
    tf = x(end-1);              % Final time of flight  
                              
    switch (cost)
        case 'Minimum time'
            switch (dynamics)
                case 'Sundman'
                    % Boundary conditions
                    P = reshape(x(1:end-2), [length(n), max(n)+1]);          % Control points
                    N = floor(x(end));                                       % The optimal number of revolutions
                    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

                    C = evaluate_state(P,B,n);                               % State evolution
                    r = sqrt(C(1,:).^2+C(3,:).^2);                           % Radial evolution

                    % Cost function
                    r = tf*trapz(tau,r);

                case 'Kepler'
                    % Cost function
                    r = tf;
            end

        case 'Minimum energy'
            % Minimize the control input
            P = reshape(x(1:end-2), [length(n), max(n)+1]);          % Control points
            N = floor(x(end));                                       % The optimal number of revolutions
            P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

            C = evaluate_state(P,B,n);                               % State evolution
            [u, ~, ~] = acceleration_control(mu, tau, C, tf, dynamics); 
            u = u/tf;
        
            % Control cost
            switch (dynamics)
                case 'Sundman'
                    r = sqrt(C(1,:).^2+C(3,:).^2);                   % Radial evolution
                    u = u./(r*trapz(tau,r));                         % Dimensional acceleration

                case 'Kepler'
            end

            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);                 % Non-dimensional acceleration norm
            
            % Cost function
            r = trapz(tau,a); 
            
        otherwise 
            error('No valid cost function was selected to be minimized');
    end       
end
