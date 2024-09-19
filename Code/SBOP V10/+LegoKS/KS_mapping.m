%% KS mapping %%
% Function to compute the mapping between the Cartesian state vector and
% the u space

% Inputs: - array x, either the Cartesian of the u-space state vector
%         - boolean direction, to determine the sense of the transformation
%         - string sundman_n, the Sundman transformation to be used in the
%           regularization
%         - scalar mu, the gravitational parameter of the system in use

% Outputs: - array S, the transformed state vector

function [S] = KS_mapping(x, direction, sundman_n, mu)

    % Sanity checks and default values 
    if ~exist("sundman_n", "var")
        sundman_n = "1";
    end

    if ~exist("mu", "var")
        mu = 1;
    end

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

        % Dimensionalize the velocity as a function of the Sundman
        % transformation
        switch sundman_n
            case "1"
                % Nothing needs to be done here
    
            case "Ecc"
                [~, alpha] = LegoKS.OscEnergy(mu, x, "Cart");                   % Energy value
                beta = sqrt(mu * alpha);                                        % Dimension factor
                S(5:8,:) = S(5:8,:) ./ beta;                                    % Dimensional velocity
    
            case "Nu"
                r = dot(S(1:4,:), S(1:4,:));                                    % Radius vector
                S(5:8,:) = S(5:8,:) .* r;                                       % Dimensional velocity
    
            otherwise
                error("The selected regularization does not have an osculatign energy equation... Aborting");
        end
    else
        % Preallocation 
        S = zeros(6, size(x,2)); 

        % Transform the position space 
        S(1:3,:) = LegoKS.KSmap(x(1:4,:), true);

        % Transformation from the U space to the Cartesian space
        r = dot(x(1:4,:), x(1:4,:));

        for i = 1:size(x,2)
            L = LegoKS.KSmatrix( x(1:4,i) );            
            aux = 2 / r(i) * L * x(5:8,i);               
            S(4:6,i) = aux(1:3);                                     
        end

        % Dimensionalizing the velocity depending on the Sundman
        % transformation
        switch sundman_n
            case "1"
                % Nothing needs to be done here
    
            case "Ecc"
                [~, alpha] = LegoKS.OscEnergy(mu, x, sundman_n);                % Energy value
                beta = sqrt(mu * alpha);                                        % Dimension factor
                S(4:6,:) = S(4:6,:) .* beta;                                    % Dimensional velocity
    
            case "Nu"
                S(4:6,:) = S(4:6,:) ./ r;                                       % Dimensional velocity
    
            otherwise
                error("The selected regularization does not have an osculatign energy equation... Aborting");
        end
    end
end