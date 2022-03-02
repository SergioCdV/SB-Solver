%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: 

% Outputs: 

function [tf] = flight_time(P, B, m, n)
    % Trajectory approximation
    P = reshape(P(1:end-m), [2, max(n)+1]);
    C = evaluate_state(P, B, n);
    
    [r,v] = extract_coordinates(C);
    v = sqrt(v(1,:).^2 + v(2,:).^2);

    % Time steps estimations along the discretization points
    dt = zeros(1,m-1);          % Preallocation for speed

    for i = 2:m
       dr = norm(r(:,i)-r(:,i-1));
       dt(i) = dr / norm(v(:,i));
    end

    % Final estimation of the time of flight
    tf = sum(dt);
end