%% Bilinear function %%
% Function to compute the bilinear operator between to vectors in U4

% Inputs: - vector u, one u vector 
%         - vector v, one u vector

% Outputs: - l, the bilinear operator between the two

function [l] = bilinear_function(u,v)
    % Compute the bilinear operator 
    l = u(1,:).*v(4,:)-u(2,:).*v(3,:)+u(3,:).*v(2,:)-u(4,:).*v(1,:);
end