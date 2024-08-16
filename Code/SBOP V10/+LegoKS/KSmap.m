%% KS map %%
% Function to compute the Hopf map between S3 and S1xS2 (and its inverse)

% Inputs: - vector u, the state variable

% Outputs: - array L, the KS operator

function [u] = KSmap(r, direction)
    % Sanity checks 
    if ~exist("direction", "var")
        direction = true;
    end

    if ( direction )
        % Pre-allocation 
        u = zeros(4, size(r,2)); 

        % Transformation from the KS space to the Cartesian space
        for i = 1:size(r,2)
            L = LegoKS.KSmatrix( r(:,i) );
            u(:,i) = L * r(:,i);
        end

        u = u(1:3,:);
        
    else
        % Pre-allocation 
        u = zeros(4, size(r,2)); 
    
        % Compute the mapping to the Hopf fibre 
        r_norm = sqrt( dot(r,r,1) ); 
        r_pos = r(1,:) >= 0; 
        r_neg = ~r_pos;

        u(4,r_pos) = zeros(1, sum(r_pos));
        u(1,r_pos) = sqrt( 0.5 * ( r(1,r_pos) + r_norm(r_pos) ) );
        u(2,r_pos) = 0.5 * ( r(2,r_pos) ./ u(1,r_pos) );
        u(3,r_pos) = 0.5 * ( r(3,r_pos) ./ u(1,r_pos) );

        u(3,r_neg) = zeros(1, sum(r_neg));
        u(2,r_neg) = sqrt( 0.5 * ( r_norm(r_neg) - r(1,r_neg) ) );
        u(1,r_neg) = 0.5 * ( r(2,r_neg) ./ u(2,r_neg) );
        u(4,r_neg) = 0.5 * ( r(3,r_neg) ./ u(2,r_neg) );
    end
end