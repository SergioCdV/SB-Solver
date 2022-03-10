%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [J] = LQR(x,B,n,m,tf)
    % Obtain the control points 
    P = reshape(x(1:end-m-1), [1 max(n)+1]);
    u = reshape(x(end-m:end-1), [1 m]); 

    % Evaluate the evolution 
    C = evaluate_state(P,B,n);

    % Cost function 
    J = (1/2)*sum((C(1,:).^2+u.^2));
end