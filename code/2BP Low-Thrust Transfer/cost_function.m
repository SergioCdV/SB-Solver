%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [r] = maximum_radius(x, B, m, n, time)
    % Maximize the orbital radius transfer
    P = reshape(x(1:end-1-2*m), [4, max(n)+1]);
    u = reshape(x(end-2*m:end-1), [2 m]);

    % Boundary conditions
    C = evaluate_state(P,B,n);

    % Maximize the radius 
    r = trapz(time*x(end), sqrt(u(1,:).^2+u(2,:).^2));
end