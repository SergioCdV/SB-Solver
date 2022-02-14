%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector s, the vector to be transformed
%         - scalar T, the dimensionalising characteristic time 
%         - scalar order, the order of the derivative to be dimensionalise
%         - vector n, the degrees of approximation to each position
%           coordinate

% Outputs: - vector r, the dimensional position vector 
%          - vector v, the dimensional velocity vector 
%          - vector gamma, the dimensional velocity vector

function [tf] = flight_time(P, B, m, T, n)
    % Trajectory approximation
    C = evaluate_state(P, B, n);
    
    % Trajectory position, velocity and acceleration coordinates
    [rho, v, ~] = extract_coordinates(C);
    v = sqrt(v(1,:).^2 + (rho(1,:).*v(2,:)).^2 + v(3,:).^2);

    % Time steps estimations along the discretization points
    dt = zeros(1,m-1);          % Preallocation for speed

    for i = 2:m
       dr = sqrt((rho(1,i)-rho(1,i-1))^2 + (rho(3,i)-rho(3,i-1))^2 + (rho(1,i)*(rho(2,i)-rho(2,i-1)))^2);
       dt(i) = dr / norm(v(:,i));
    end

    % Final estimation of the time of flight
    tf = sum(dt)*T;
end