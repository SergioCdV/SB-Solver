%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial trajectory guess given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - scalar n, the order of the approximating polynomail function
%         - vector x0, the initial 6 by 1 state vector (position and velocity)
%         - vector xf, the final 6 by 1 state vector (position and velocity)
%         - array P0, the initial array of control points 
%         - cell array B, the polynomial basis to be used
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the updated boundary conditions control points, of dimensions size(x,1) x n+1 

function [P] = boundary_conditions(n, x0, xf, P0, B, basis)
    % Sanity check 
    if (min(n) < 2)
        error('Cauchy boundary conditions cannot be imposed');
    else
        P = P0;         % Initialization
    end

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = x0;
            for i = 1:length(n)
                P(i,n(i)+1) = xf(i);
            end

        case 'Chebyshev'
            % Dimensionality check 
            if (size(x0,2) ~= 1)
                x0 = x0.';
            end
        
            if (size(xf,2) ~= 1)
                xf = xf.';
            end
        
            for i = 1:size(P,1)
                % Symbolic regression 
                l = n(i)-2;
                X0(1) = x0(i)-sum(P0(i,3:n(i)-1).*(-1).^(2:l),2); 
                Xf(1) = xf(i)-sum(P0(i,3:n(i)-1),2);

                A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); ...
                     1 1 1 1];

                sol = A\[X0(1); Xf(1)];
                P(i,[1 2 n(i) n(i)+1]) = sol.';
            end

        case 'Legendre'
            % Dimensionality check 
            if (size(x0,2) ~= 1)
                x0 = x0.';
            end
        
            if (size(xf,2) ~= 1)
                xf = xf.';
            end

            for i = 1:size(P,1)
                % Symbolic regression 
                l = n(i)-2;
                X0(1) = x0(i)-sum(P0(i,3:n(i)-1).*(-1).^(2:l),2); 
                Xf(1) = xf(i)-sum(P0(i,3:n(i)-1),2);

                A = [1 -1 (-1)^(n(i)-1) (-1)^n(i); ...
                     1 1 1 1];

                sol = A\[X0(1); Xf(1)];
                P(i,[1 2 n(i) n(i)+1]) = sol.';
            end

        otherwise
            % Compute the partial state evolution 
            N = size(P,1);                   % Number of state variables
            C = zeros(2*N,2);                % Preallocation for speed
            for i = 1:length(n)
                C(i,:) = P0(i,3:n(i)-1)*B{i}(3:n(i)-1,[1 end]);
                C(length(n)+i,:) = P0(i,3:n(i)-1)*B{i}(n(i)+4:2*n(i),[1 end]);
            end

            % Assemble the linear system 
            C = [x0.'-C(1:N,1) xf.'-C(1:N,end)];
            
            % Control points for an orthogonal Bézier curve
            for i = 1:length(n)
                index = [1 n(i)+1];
                A = B{i}(index,[1,end]);
                P(i,[1 2 n(i) n(i)+1]) = C(i,:)*A^(-1);
            end
    end
end