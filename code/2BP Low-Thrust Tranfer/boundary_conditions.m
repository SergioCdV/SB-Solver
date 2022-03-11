%% Project: 
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial Bézier points given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - vector x0, the initial 6 by 1 state vector (heliocentric position and velocity)
%         - vector xf, the final 6 by 1 state vector (heliocentric position and velocity)
%         - scalar tfapp, the approximated time of flight initial guess
%         - scalar n, the order of the approximating Bézier curve
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the boundary conditions control points, of dimensions 3 x n+1 

function [P] = boundary_conditions(mu, tfapp, n, x0, xf, basis)
    % Constants 
    P = zeros(2,2);            % Preallocation of the boundary control points

    % Correct dimensionalising
    x0(4) = x0(4)/x0(1);
    xf(4) = xf(4)/xf(1);

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = x0(1:2);            
            P(:,2) = x0(1:2) + x0(3:4)*tfapp./n;          
            P(:,3) = xf(1:2) - xf(3:4)*tfapp./n;
            P(:,4) = xf(1:2);

        case 'Orthogonal Bernstein'
            % Control points for an orthogonal Bézier curve
            P(:,1) = x0(1:2);            
            P(:,2) = x0(1:2) + x0(3:4)*tfapp./n;          
            P(:,3) = xf(1:2) - xf(3:4)*tfapp./n;
            P(:,4) = xf(1:2);

        otherwise 
            error('No valid collocation polynomial basis has been selected');
    end
end