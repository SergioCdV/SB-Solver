%% Project: SBOPT %%
% Date: 05/05/2023

%% Grid %% 
% Function to compute a desired polynomail collocation grid

% Inputs: - scalar m, the number of points in the grid

% Outputs: - object Grid, defining the collocation grid and associated quadrature rule

function [Grid] = gridding(obj, m)
    
    if (~exist('m', 'var'))
        m = obj.NumNodes;
    end

    % Final sampling distribution setup
    switch (obj.Grid)
        
        case 'Chebyshev'
            Grid = CollocationMesh.ChebyshevGrid(m);

        case 'Legendre'
            Grid = CollocationMesh.LegendreGrid(m);

        case 'Linear'
            Grid = CollocationMesh.LinearGrid(m);

        case 'Normal'
            Grid = CollocationMesh.NormalGrid(m); 

        case 'Random'
            Grid = CollocationMesh.RandomGrid(m); 

        case 'Newton-Cotes'
            Grid = CollocationMesh.NewtonCotesGrid(m); 

        case 'Trapezoidal'
            Grid = CollocationMesh.TrapezoidalGrid(m);

        case 'Bernstein'
            Grid = CollocationMesh.BezierGrid(m);

        case 'Orthogonal Bernstein'
            Grid = CollocationMesh.OBezierGrid(m);

        otherwise
            error('No valid quadrature was selected');
    end
end