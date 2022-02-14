%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to dimensionalise a given set of coodinates

% Inputs: - vector s, the vector to be transformed

% Outputs: - vector r, the dimensional position vector 
%          - vector v, the dimensional velocity vector 
%          - vector gamma, the dimensional velocity vector

function  [r, v, gamma] = extract_coordinates(s)
    % Dimensional position coordinates 
      r = s(1:3,:);

    % Dimensional velocity coordinates 
      v = s(4:6,:);

    % Dimensional acceleration coordinates
      gamma = s(7:9,:);
end