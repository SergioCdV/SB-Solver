%% Project: 
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial Bézier points given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - scalar n, the order of the approximating Bézier curve
%         - vector t, the time measurements vector
%         - array r, the position measurements array
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the boundary conditions control points, of dimensions 3 x n+1 

function [P] = boundary_conditions(n, t, r, basis)
    % Constants 
    P = zeros(size(r,1),2);            % Preallocation of the boundary control points

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = r(:,1);
            P(:,2) = r(:,end);

        case 'Orthogonal Bernstein'
            % Assemble the linear system 
            P(:,1) = r(:,1);
            P(:,2) = r(:,end);

        otherwise 
            error('No valid collocation polynomial basis has been selected');
    end
end