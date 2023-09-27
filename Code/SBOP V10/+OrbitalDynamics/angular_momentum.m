%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 18/04/22

%% Angular momentum %%
% Compute the angular momentum vector evolution given a cylindrical state evolution C

% Inputs: - array C, the cylindrical state evolution
%         - string method, the parameter distribution to be used

% Outputs: - vector h, the angular momentum norm evolution

function [h] = angular_momentum(C)
    % Preallocation 
    h = zeros(1,size(C,2)); 

    % Compute the angular momentum
    for i = 1:size(C,2)
        s = OrbitalDynamics.cylindrical2cartesian(C(1:6,i), true);
        h(i) = norm(cross(s(1:3),s(4:6)));
    end
end