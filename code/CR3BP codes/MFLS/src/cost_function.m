%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 18/08/22
% File: mfls_optimization.m 
% Issue: 0 
% Validated: 18/08/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - string cost, indicating the cost function to be minimized 
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(wp, wv, kap, phi, psi, cost, initial, final, n, tau, x, B, basis)
    % Optimization variables 
    tf = x(end);                                      % The final time of flight

    switch (cost)
        case 'Minimum time'
            r = tf;
            
        case 'Minimum energy'
            % Minimize the control input
            P = reshape(x(1:end-1), [length(n), max(n)+1]);     % Control points
        
            % Boundary conditions
            P = boundary_conditions(tf, n, initial, final, P, B, basis);
        
            % State evolution
            C = evaluate_state(P,B,n);
        
            % Control input
            u = acceleration_control(wp, wv, kap, phi, psi, C, tf, tau);        
        
            % Control cost
            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Dimensional acceleration

            % Cost function
            r = trapz(tf*tau,a);
            
        otherwise 
            error('No valid cost function was selected');
    end
end