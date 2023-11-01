classdef Solver
    properties 
        % Numerical solver configuration
        PolOrder;       % Polynomial order for each state dimension
        NumNodes;       % Number of nodes in the independent variable grid
        Basis;          % Polynomial basis 
        Grid;           % Define the independent variable grid to be used

        % Initial guess 
        InitialGuessFlag = false;   % Flag to indicate an initial guess is supplied
        P0;                         % Parameters initial guess

        % Optimization configuration 
        maxIter = 1e4;
        maxFunctionEvaluations = 1e6;
    end

    methods 
        % Solver definition 
        function [obj] = Solver(myBasis, myPolOrder, myGrid, myNumNodes)
            % Solver definition
            if (exist('myBasis', 'var'))
                obj.Basis = myBasis;
            else
                obj.Basis = 'Legendre';
            end
            
            if (exist('myPolOrder', 'var'))
                obj.PolOrder = myPolOrder;
            else
                obj.PolOrder = 10;
            end
            
            if (exist('myGrid', 'var'))
                obj.Grid = myGrid;
            else
                obj.Grid = 'Legendre';
            end
    
            if (exist('myNumNodes', 'var'))
                obj.NumNodes = myNumNodes;
                else
                obj.NumNodes = ceil((obj.PolOrder - 1)/2);
            end

            % Checks and environment settings
            obj = obj.Check();
            obj.set_graphics();
        end

        % Solve
        [C, cost, u, t0, tf, t, exitflag, output, P] = solve(obj, Problem);
        [B, C, tau] = state_basis(obj, L, n, basis, tau);
        [Grid] = gridding(obj,m);
        [C] = evaluate_state(obj, P, B, n, L);
    end

    methods (Access = private)
        % Check function 
        function [obj] = Check(obj)
            % Check the numerical solver 
            if (obj.NumNodes < 2 * max(obj.PolOrder)+1)
                warning('Quadrature may not be accurate. Consider increasing the number of nodes in the grid.');
            end

            switch (obj.Basis)
                case 'Bernstein'
                    switch (obj.Grid)
                        case 'Chebyshev'
                            error('Selected grid distribution is incompatible with Bernstein polynomials');
                        case 'Legendre'
                            error('Selected grid distribution is incompatible with Bernstein polynomials');
                        otherwise
                    end
 
                case 'Orthogonal Bernstein'
                    switch (obj.Grid)
                        case 'Chebyshev'
                            error('Selected grid distribution is incompatible with orthogonal Bernstein polynomials.');
                        case 'Legendre'
                            error('Selected grid distribution is incompatible with orthogonal Bernstein polynomials.');
                        otherwise
                    end

                case 'Chebyshev'
                    switch (obj.Grid)
                        case 'Chebyshev'
                        otherwise
                            warning('Consider selecting a Clenshaw-Curtis independent grid for maximum accuracy.');
                    end

                case 'Legendre'
                    switch (obj.Grid)
                        case 'Legendre'
                        otherwise
                            warning('Consider selecting a Lagrange-Gauss-Lobatto independent grid for maximum accuracy.');
                    end

                otherwise
                   error('The selected polynomial support is not currently supported.')
            end
        end
        
        [beta, t0, tf, P, C] = initial_approximation(obj, Problem, B, Grid);
        [P, C] = initial_fitting(obj, Problem, Grid, s);
        [P_lb, P_ub] = opt_bounds(obj, Problem, n, B);
        [P] = boundary_conditions(obj, Problem, beta, t0, tf, tau, B, basis, n, P0);
        [c, ceq] = constraints(obj, Problem, B, CB, Grid, x);
        [r] = cost_function(obj, Problem, B, Grid, x);
    end

    methods (Static, Access = private)
        % Set graphics
        function set_graphics()
            % Set graphical properties
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
            set(groot, 'defaultAxesFontSize', 11); 
            set(groot, 'defaultAxesGridAlpha', 0.3); 
            set(groot, 'defaultAxesLineWidth', 0.75);
            set(groot, 'defaultAxesXMinorTick', 'on');
            set(groot, 'defaultAxesYMinorTick', 'on');
            set(groot, 'defaultFigureRenderer', 'painters');
            set(groot, 'defaultLegendBox', 'off');
            set(groot, 'defaultLegendInterpreter', 'latex');
            set(groot, 'defaultLegendLocation', 'best');
            set(groot, 'defaultLineLineWidth', 1); 
            set(groot, 'defaultLineMarkerSize', 3);
            set(groot, 'defaultTextInterpreter','latex');
        end

        display_results(exitflag, cost, output);
    end
end