%% Project: 
% Date: 30/01/2022

%% Cylindrical to cartesian %% 
% Function to transform state vector between cartesian to cyclindral
% coordinates or viceversa

% Inputs: - vector s, the state vector to be transformed 
%         - boolean direction, 0 for the cylindrical to cartesian
%           transformation and 1 for the viceversa case

% Outputs: - vector S, the transformed state vector

function [S] = cylindrical2cartesian(s, direction)
    % Sanity check on the s dimensions 
    if (size(s,1) == 3)
        lastwarn('State vector has only 3 dimensions')
        s = [s; zeros(3,size(s,2))];
    end

    % Switch directions 
    if (direction)
        % Cylindrical position coordinates
        rho = s(1,:);
        theta = s(2,:);
        z = s(3,:);

        % Cylindrical velocity coordinates
        drho = s(4,:);
        dtheta = s(5,:);
        dz = s(6,:);

        % Cartesian position coordinates
        S(1,:) = rho.*cos(theta);
        S(2,:) = rho.*sin(theta);
        S(3,:) = z;

        % Cartesian velocity coordinates 
        S(4,:) = drho-rho.*sin(theta).*dtheta;
        S(5,:) = drho+rho.*cos(theta).*dtheta;
        S(6,:) = dz;
    
    else
        % Cartesian position coordinates
        x = s(1,:);
        y = s(2,:);
        z = s(3,:);

        % Cartesian velocity coordinates
        dx = s(4,:);
        dy = s(5,:);
        dz = s(6,:);

        % Cartesian position coordinates
        S(1,:) = sqrt(x.^2+y.^2);
        S(2,:) = atan2(y,x);
        S(3,:) = z;

        % Cartesian velocity coordinates 
        S(4,:) = (x.*dx+y.*dy)./S(1,:);
        S(5,:) = (x.*dy-y.*dx)./S(1,:).^2;
        S(6,:) = dz;
    end  
end