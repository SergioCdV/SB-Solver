%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 17/11/22

%% u mapping %%
% Function to compute the mapping from the position space to the u space

% Inputs: - vector u, the state variable

% Outputs: - array L, the KS operator

function [u] = u_mapping(r)
    % Pre-allocation 
    u = zeros(4,1); 

    % Compute the mapping to the Hopf fibre 
    if (r(1) >= 0)
        u(4) = 0;
        u(1) = ((r(1)+norm(r))/2)^(1/2);
        u(2) = (1/2)*(r(2)/u(1));
        u(3) = (1/2)*(r(3)/u(1));
    else
        u(3) = 0;
        u(2) = ((norm(r)-r(1))/2)^(1/2);
        u(1) = (1/2)*(r(2)/u(2));
        u(4) = (1/2)*(r(3)/u(2));
    end

        u(4) = 0;
        u(1) = ((r(1)+norm(r))/2)^(1/2);
        u(2) = (1/2)*(r(2)/u(1));
        u(3) = (1/2)*(r(3)/u(1));
end