%% Project: 
% Date: 30/01/22

%% Acceleration %%
% Function to compute the acceleration vector norm from cylindrical coordinates

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar tfapp, the approximated time of flight of the mission
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector n, ranging the polynomial basis cardinal

% Outputs: - scalar a, the acceleration vector norm

function [a] = acceleration(mu, r0, tfapp, P, B, n)
    % initialize cell array
    C = cell(1,3);
    
    % Extract the spacecraft coodinates evaluating the Bézier curve approximation
    for i = 1:3 
        for j = 1:3 
            C{i}(:,j) = squeeze(B{i}(j,:,:))*P(i,1:(n(i)+1))';
        end
    end
        
    [rho, theta, z] = extract_coordinates(C, r0, tfapp);
    
    % Heliocentric position vector
    r = sqrt(rho.o.^2 + z.o.^2);
    
    % Equations of motion
    arho = rho.DD - rho.o.*theta.D.^2 + mu.*rho.o./r.^3;
    atheta = rho.o.*theta.DD + 2.*rho.D.*theta.D;
    az = z.DD + mu.*z.o./r.^3;
    
    % Magnitude of the acceleration
    a = sqrt(arho.^2 + atheta.^2 + az.^2);
end