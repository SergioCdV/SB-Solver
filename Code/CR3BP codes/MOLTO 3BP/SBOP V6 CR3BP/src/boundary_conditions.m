%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Boundary conditions %% 
% Function to compute the initial trajectory guess given the boundary
% conditions of the spacecraft in non-dimensional units

% Inputs: - scalar tfapp, the approximated time of flight initial guess
%         - scalar n, the order of the approximating polynomail function
%         - vector x0, the initial 6 by 1 state vector (position and velocity)
%         - vector xf, the final 6 by 1 state vector (position and velocity)
%         - scalar thetaf, the final optimal insertion phase
%         - array P0, the initial array of control points 
%         - cell array B, the polynomial basis to be used
%         - string basis, to select the polynomial basis to be used in the
%           approximation

% Outputs: - array P, the updated boundary conditions control points, of dimensions size(x,1) x n+1 

function [P] = boundary_conditions(tfapp, n, x0, xf, thetaf, P0, B, basis)
    % Sanity check 
    if (min(n) < 3)
        error('Cauchy boundary conditions cannot be imposed');
    else
        P = P0;         % Initialization
    end

    % Add the revolutions to the final angle
    % xf(2) = thetaf;

    % Dimensionalizing of the velocity 
    x0(4:6) = tfapp*x0(4:6);
    xf(4:6) = tfapp*xf(4:6);

    % Switch the polynomial basis to be used
    switch (basis)
        case 'Bernstein'                
            % Control points for a nonorthogonal Bézier curve
            P(:,1) = x0(1:3);
            P(:,2) = x0(1:3)+x0(4:6)./n;
            for i = 1:length(n)
                P(i,n(i)) = xf(i)-xf(length(xf)/2+i)/n(i);
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
        
            % Symbolic regression 
            l = size(P0,2)-3;
            x0(1:3) = x0(1:3)-sum(P0(:,3:end-2).*(-1).^(2:l),2); 
            xf(1:3) = xf(1:3)-sum(P0(:,3:end-2),2);
            x0(4:6) = x0(4:6)-sum((-1).^(1:(l-1)).*P0(:,3:end-2).*(2:l).^2,2);
            xf(4:6) = xf(4:6)-sum(P0(:,3:end-2).*(2:l).^2,2);

            A = [1 -1 (-1)^(size(P0,2)-1) (-1)^size(P0,2); ...
                 1 1 1 1; ... 
                 0 1 (-1)^(size(P0,2)-3)*(size(P0,2)-2)^2 (-1)^(size(P0,2)-2)*(size(P0,2)-1)^2; ...
                 0 1 (size(P0,2)-2)^2 (size(P0,2)-1)^2];

            for i = 1:size(P,1)
                sol = A\[x0(i); xf(i); x0(3+i); xf(3+i)];
                P(i,[1 2 end-1 end]) = sol.';
            end
            
        case 'Legendre'
            % Dimensionality check 
            if (size(x0,2) ~= 1)
                x0 = x0.';
            end
        
            if (size(xf,2) ~= 1)
                xf = xf.';
            end

            % Symbolic regression 
            l = size(P0,2)-3;
            x0(1:3) = x0(1:3)-sum(P0(:,3:end-2).*(-1).^(2:l),2); 
            xf(1:3) = xf(1:3)-sum(P0(:,3:end-2),2);
            x0(4:6) = x0(4:6)-sum((-1).^(1:(l-1)).*P0(:,3:end-2).*(2:l).*((3:l+1)/2),2);
            xf(4:6) = xf(4:6)-sum(P0(:,3:end-2).*(2:l).*((3:l+1)/2),2);

            A = [1 -1 (-1)^(size(P0,2)-1) (-1)^size(P0,2); ...
                 1 1 1 1; ... 
                 0 1 (-1)^(size(P0,2)-3)*(size(P0,2)-2)*(size(P0,2)-1)/2 (-1)^(size(P0,2)-2)*(size(P0,2)-1)*(size(P0,2))/2; ...
                 0 1 (size(P0,2)-2)*(size(P0,2)-1)/2 (size(P0,2)-1)*(size(P0,2))/2];

            for i = 1:size(P,1)
                sol = A\[x0(i); xf(i); x0(3+i); xf(3+i)];
                P(i,[1 2 end-1 end]) = sol.';
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
            C = [x0.'-C(:,1) xf.'-C(:,end)];
            C = [C(1:size(P,1),:) C(size(P,1)+1:2*size(P,1),:)];
            
            % Control points for an orthogonal Bézier curve
            for i = 1:length(n)
                index = [1 2 n(i) n(i)+1 n(i)+2 n(i)+3 2*n(i)+1 2*(n(i)+1)];
                A = B{i}(index,[1,end]);
                A = [A(1:4,:) A(5:8,:)];
                P(i,[1 2 n(i) n(i)+1]) = C(i,:)*A^(-1);
            end
    end
end