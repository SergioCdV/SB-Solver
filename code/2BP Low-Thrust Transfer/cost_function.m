%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [r] = cost_function(x, B, m, n, time)
    % Maximize the orbital radius transfer
    P = reshape(x(1:end-1-3*m), [6, max(n)+1]);
    u = reshape(x(end-3*m:end-1), [3 m]);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    % Maximize the radius 
    r = x(end)*trapz(time, u(1,:).^2+u(2,:).^2+u(3,:).^2);
end