%% Project: 
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial Bézier points given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - vector x0, the initial 6 by 1 state vector (heliocentric position and velocity)
%         - vector xf, the final 6 by 1 state vector (heliocentric position and velocity)
%         - scalar n, the order of the approximating Bézier curve
%         - scalar r0, the characteristic or dimensionalising distance 
%         - scalar T, the characteristic or dimensionalising time
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the boundary conditions control points, of dimensions 3 x n+1 

function [P] = boundary_conditions(x0, xf, n, r0, T, basis)
    % Constants 
    adimv = [r0; 1; r0];        % Dimensionalising units vector
    P = zeros(3, 4);            % Preallocation of the boundary control points

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = x0(1:3)./adimv;            
            P(:,2) = x0(1:3)./adimv + x0(4:6)*T./n./adimv;          
            P(:,3) = xf(1:3)./adimv - xf(4:6)*T./n./adimv;
            P(:,4) = xf(1:3)./adimv;

        case 'Orthogonal Bernstein'
            % Control points for an orthogonal Bézier curve
            P(:,1) = x0(1:3)./adimv;            
            P(:,2) = x0(1:3)./adimv + x0(4:6)*T./n./adimv;          
            P(:,3) = xf(1:3)./adimv - xf(4:6)*T./n./adimv;
            P(:,4) = xf(1:3)./adimv;

        otherwise 
            error('No valid collocation polynomial basis has been selected');
    end
end