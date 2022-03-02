%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to dimensionalise a given set of coodinates

% Inputs: - vector s, the vector to be transformed

% Outputs: - vector r, the dimensional position vector 
%          - vector v, the dimensional velocity vector

function  [r, v] = extract_coordinates(s)
    % Preallocation 
    r = zeros(2,size(s,2));         % Position vector 
    v = zeros(3,size(s,2));         % Velocity vector 

    % Dimensional position coordinates 
    r(1,:) = s(1,:);
    r(2,:) = s(2,:);

    % Dimensional velocity coordinates 
    v(1,:) = s(3,:);
    v(2,:) = s(4,:);
end