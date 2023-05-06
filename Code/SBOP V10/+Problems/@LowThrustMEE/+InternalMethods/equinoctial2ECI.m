%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Classical orbital elements to equinoctial %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - scalar mu, the gravitational parameter of the system
%         - vector s, the orbital state vector to be transformed 
%         - boolean direction, 0 for the COE to equinoctional
%           transformation and 1 for the viceversa case

% Outputs: - vector S, the transformed orbital state vector

function [S] = equinoctial2ECI(obj, mu, s, direction)
    % Sanity check on the s dimensions 
    if (size(s,1) ~= 6)
        lastwarn('Input orbital element set is not valid')
    end

    % Switch directions 
    if (direction)
        % Auxiliary variables 
        alpha = s(4,:).^2-s(5,:).^2;
        m = 1+s(4,:).^2+s(5,:).^2;
        w = 1+s(2,:).*cos(s(6,:))+s(3,:).*sin(s(6,:));
        r = s(1,:)./w;

        % MEE to ECI conversion 
        S(1,:) = r./m.*(cos(s(6,:))+alpha.*cos(s(6,:))+2*s(4,:).*s(5,:).*sin(s(6,:)));
        S(2,:) = r./m.*(sin(s(6,:))-alpha.*sin(s(6,:))+2*s(4,:).*s(5,:).*cos(s(6,:)));
        S(3,:) = 2*r./m.*(s(4,:).*sin(s(6,:))-s(5,:).*cos(s(6,:)));
        S(4,:) = -1./m.*sqrt(mu./s(1,:)).*(sin(s(6,:))+alpha.*sin(s(6,:))-2*s(4,:).*s(5,:).*cos(s(6,:))+s(3,:)-2*s(2,:).*s(4,:).*s(5,:)+alpha.*s(3,:));
        S(5,:) = -1./m.*sqrt(mu./s(1,:)).*(-cos(s(6,:))+alpha.*cos(s(6,:))+2*s(4,:).*s(5,:).*sin(s(6,:))-s(2,:)+2*s(3,:).*s(4,:).*s(5,:)+alpha.*s(2,:));
        S(6,:) = 2./m.*sqrt(mu./s(1,:)).*(s(4,:).*cos(s(6,:))+s(5,:).*sin(s(6,:))+s(4,:).*s(2,:)+s(5,:).*s(3,:));
    
    else
        % Preallocation 
        S = zeros(size(s));

        % ECI to MEE transformation
        for i = 1:size(s,2)
            coe = state2coe(mu, s(:,i), 'Inertial');
            S(:,i) = coe2equinoctial(coe,true);
        end
    end  
end