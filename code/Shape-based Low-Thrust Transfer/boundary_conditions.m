%% Project: 
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial Bézier points given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - scalar tfapp, the approximated time of flight initial guess
%         - scalar n, the order of the approximating Bézier curve
%         - vector x0, the initial 6 by 1 state vector (heliocentric position and velocity)
%         - vector xf, the final 6 by 1 state vector (heliocentric position and velocity)
%         - scalar N, the number of revolutions 
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the boundary conditions control points, of dimensions 3 x n+1 

function [P] = boundary_conditions(tfapp, n, x0, xf, N, basis)
    % Constants 
    P = zeros(length(x0)/2,4);            % Preallocation of the boundary control points

    % Add the revolutions to the final angle
    xf(2) = xf(2)+2*pi*N;

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = x0(1:3);
            P(:,2) = x0(1:3)+tfapp*x0(4:6)./n;
            P(:,3) = xf(1:3)-tfapp*xf(4:6)./n;
            P(:,4) = xf(1:3);

        case 'Orthogonal Bernstein'
            % Assemble the linear system 
            b = [x0 xf].';
            b = [b([1 7 4 10]); b([2 8 5 11]); b([3 9 6 12])];
            A = [OB_basis(n(i), [0 1]).'; OB_derivative(n(i), [0 1], 1).']; 
            A = kron(eye(3),A);

            % Control points for an orthogonal Bézier curve
            P = A\b; 
            P = reshape(P, [n 3]).';

        otherwise 
            error('No valid collocation polynomial basis has been selected');
    end
end