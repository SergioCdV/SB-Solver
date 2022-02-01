%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector s, the vector to be transformed
%         - scalar r0, the dimensionalising characteristic distance 
%         - scalar tf, the estimated time of flight 
%         - scalar order, the order of the derivative to be dimensionalise

% Outputs: - vector r, the dimensional position vector 
%          - vector v, the dimensional velocity vector 
%          - vector gamma, the dimensional velocity vector

function [tf] = flight_time(P, B, m, tf, r0)
    % Trajectory approximation
    C = zeros(9,size(B,2));     % Preallocation for speed
    k = size(B,1)/3;            % Number of control points
    for i = 1:3
        C(1+3*(i-1):3*i,:) = P*B(1+k*(i-1):k*i,:);
    end
    
    % Trajectory position, velocity and acceleration coordinates
    [rho, v, ~] = extract_coordinates(C, r0, tf);
    v = sqrt(v(1,:).^2 + (rho(1,:).*v(2,:)).^2 + v(3,:).^2);

    % Time steps estimations along the discretization points
    dt = zeros(1,m-1);          % Preallocation for speed

    for i = 2:m
       dr = sqrt((rho(1,i)-rho(1,i-1))^2 + (rho(3,i)-rho(3,i-1))^2 + (rho(1,i)*(rho(2,i)-rho(2,i-1)))^2);
       dt(i) = dr / norm(v(:,i));
    end

    % Final estimation of the time of flight
    tf = sum(dt);
end