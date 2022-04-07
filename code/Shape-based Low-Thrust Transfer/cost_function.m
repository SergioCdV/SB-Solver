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
    P = reshape(x(1:end-1), [length(n), max(n)+1]);
    tf = x(end);

    % Boundary conditions
    P(:,1) = initial(1:3);
    P(:,2) = initial(1:3)+tf*initial(4:6)./n;
    P(:,end-1) = final(1:3)-tf*final(4:6)./n;
    P(:,end) = final(1:3);

    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(mu,C,tf); 
    u = u/tf^2;

    % Minimize the control input
    r = trapz(time*tf, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2));
end