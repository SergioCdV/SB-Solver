%% Project:  
% Sergio Cuevas del Valle
% Date: 20/01/20
% File: OB_nodes
% Issue: 0 
% Validated: 

%% Orthogonoal Bézier nodes %%
% This scripts provides the function to compute the orthogonal Bézier nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest
%         - scalar a, the initial domain boundary 
%         - scalar b, the final domain boundary 

% Output: - vector y, the Bézier nodes of interest

function [y] = OB_nodes(a,b,N)
    % Polynomial order 
    n = N+1; 

    % Initial guess for the nodes
    y = LG_nodes(a,b,N);
    y = linspace(0,1,N);

    % Main loop
    for i = 1:N
        % Newton method setup 
        GoOn = true; 
        iter = 0; 
        maxIter = 100; 
        tol = 1e-3; 

        % Initial guess 
        roots = y(i);

        % Newton method 
        while (GoOn && iter < maxIter)
            % Compute the Bézier polynomial and its derivative 
            f = OB_basis(i, roots);   
            df = OB_derivative(i, roots, 1);

            % Newton step 
            ds = -f(end)/df(end);
            roots = roots + ds;

            % Convergence analysis 
            if (abs(ds) < tol)
                GoOn = false;
            else
                iter = iter+1; 
            end
        end

        % Save the results 
        y(i) = roots;
    end
end