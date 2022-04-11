%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [r] = cost_function(initial, final, mu, T, x, B, m, n, time)
    % Minimize the control input
    P = reshape(x(1:end-2), [length(n), max(n)+1]);
    tf = x(end-1);
    N = floor(x(end));

    % Boundary conditions
    P(:,[1 2 end-1 end]) = boundary_conditions(mu, tf, n, initial, final, N, 'Orthogonal Bernstein');

    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(mu,C,tf); 
    a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);

    % Minimize the control input
    r = trapz(time, a/tf);
end