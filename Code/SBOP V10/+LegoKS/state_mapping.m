%% State mapping %%
% Function to compute the mapping between the Cartesian state vector and
% the u space

% Inputs: - array x, either the Cartesian of the u-space state vector
%         - boolean direction, to determine the sense of the transformation

% Outputs: - array S, the transformed state vector

function [S] = state_mapping(x, direction)

    if ( direction )
        % Preallocation 
        S = zeros(8, size(x,2)); 

        % Transform the position space 
        S(1:4,:) = LegoKS.KSmap(x(1:3,:), false);

        % Transformation from the Cartesian space to the U space
        for i = 1:size(x,2)
            L = LegoKS.KSmatrix( S(1:4,i) );                          
            S(5:8,i) = 0.5 * L.' * [x(4:6,i); 0];                           
        end
    else
        % Preallocation 
        S = zeros(6, size(x,2)); 

        % Transform the position space 
        S(1:3,:) = LegoKS.KSmap(x(1:4,:), true);

        % Transformation from the U space to the Cartesian space
        for i = 1:size(x,2)
            L = LegoKS.KSmatrix( x(1:4,i) );            
            aux = 2 / dot(x(1:4,i), x(1:4,i)) * L * x(5:8,i);               
            S(4:6,i) = aux(1:3);                                     
        end
    end
end