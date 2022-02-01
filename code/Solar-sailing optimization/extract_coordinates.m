%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to dimensionalise a given set of coodinates

% Inputs: - vector s, the vector to be transformed
%         - scalar r0, the dimensionalising characteristic distance 
%         - scalar tf, the estimated time of flight 
%         - scalar order, the order of the derivative to be dimensionalise

% Outputs: - vector r, the dimensional position vector 
%          - vector v, the dimensional velocity vector 
%          - vector gamma, the dimensional velocity vector

function  [r, v, gamma] = extract_coordinates(s, r0, tf)
    % Preallocation 
    r = zeros(3,size(s,2));         % Position vector 
    v = zeros(3,size(s,2));         % Velocity vector 
    gamma = zeros(3,size(s,2));     % Acceleration vector 

    % Non-dimensional position coordinates 
    r(1,:) = s(1,:)*r0;
    r(2,:) = s(2,:);
    r(3,:) = s(3,:)*r0;

    % Non-dimensional velocity coordinates 
    v(1,:) = s(4,:)*r0/tf;
    v(2,:) = s(5,:)/tf;
    v(3,:) = s(6,:)*r0/tf;

    % Non-dimensional acceleration coordinates
    gamma(1,:) = s(7,:)*r0/tf^2;
    gamma(2,:) = s(8,:)/tf^2;
    gamma(3,:) = s(9,:)*r0/tf^2;
end