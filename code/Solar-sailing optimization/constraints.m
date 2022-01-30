%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar tfapp, the approximated time of flight of the mission
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector n, ranging the polynomial basis cardinal
%         - scalar amax, the maximum allowed acceleration magnitude

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, r0, tfapp, P, P0, B, n, amax)
    % Non-linear inequality on the acceleration magnitude (a < a_max)
    amag = acceleration(mu, r0, tfapp, P, B, n);
    c = amag - amax;

    % Boundary conditions
    ceq = zeros(3,2);

    for k = 1:3
        L = n(k)+1;
        P0aux = P0(k,1:L);
        Paux = P(k,1:L);
        deltaP = Paux - P0aux;
    
        ceq(k,:) = [deltaP(1),  deltaP(L)];
    end
end