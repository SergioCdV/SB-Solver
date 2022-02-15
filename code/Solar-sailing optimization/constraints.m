%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar tfapp, the approximated time of flight of the mission
%         - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - scalar amax, the maximum allowed acceleration magnitude

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, r0, n, x, B, amax)
    % Extract the optimization variables
    P = reshape(x(1:end-1), [3, max(n)+1]);
    tf = x(end);

    % Non-linear inequality on the acceleration magnitude (a < a_max)
    amag = acceleration(mu, r0, tf, P, B, n);
    c = (amag - amax*ones(size(amag)));
    c = c.';
    c = zeros(size(c));

    % Boundary conditions
    C = evaluate_state(P,B,n);
    ceq = [C(1:3,1)-initial(1:3)./[r0;1;r0]; ...
           C(1:3,end)-final(1:3)./[r0;1;r0]; ...
           C(4:6,1)-initial(4:6)./[r0;1;r0]*tf; ...
           C(4:6,end)-final(4:6)./[r0;1;r0]*tf];
end