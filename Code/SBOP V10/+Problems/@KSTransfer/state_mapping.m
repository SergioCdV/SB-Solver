%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 17/11/22

%% State mapping %%
% Function to compute the mapping between the Cartesian state vector and
% the u space

% Inputs: - array x, either the Cartesian of the u-space state vector
%         - boolean direction, to determine the sense of the transformation

% Outputs: - array S, the transformed state vector

function [S] = state_mapping(x, direction)
    % Compute the mapping
    if (direction)
        % Preallocation 
        S = zeros(8, size(x,2)); 

        % Transformation from the Cartesian space to the U space
        for i = 1:size(x,2)
            S(1:4,i) = Problems.KSTransfer.u_mapping(x(1:3,:));           % Position space transformation
            L = Problems.KSTransfer.KS_matrix(S(1:4,i));                  % KS matrix
            S(5:8,i) = (1/2)*L.'*[x(4:6,:); zeros(1,size(x,2))];          % Velocity space transformation
        end
    else
        % Preallocation 
        S = zeros(6, size(x,2));

        % Transformation from the U space to the Cartesian space
        for i = 1:size(x,2)
            L = Problems.KSTransfer.KS_matrix(x(1:4,i));             % KS matrix
            aux = L*x(1:4,i);                                        % Position space transformation
            S(1:3,i) = aux(1:3);                                     % Position space transformation
            aux = 2/dot(S(1:4,i),S(1:4,i))*L*x(5:8,:);               % Velocity space transformation
            S(4:6,i) = aux(1:3);                                     % Velocity space transformation
        end
    end
end