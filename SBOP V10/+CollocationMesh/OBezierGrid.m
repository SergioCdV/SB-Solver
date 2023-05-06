%% Project: SBOPT %%
% Date: 05/05/23

%% Bezier quadrature object %%
% Definition of the Bezier quadrature rule and collocation grid as an
% object

classdef BezierGrid < CollocationMesh.AbstractGrid
    
    methods 
        % Constructor 
        function [obj] = BezierGrid(m)
            if (exist('m', 'var'))
                obj.NumberPoints = m;
            else
                obj.NumberPoints = 100;
            end
 
            % Generation of the grid and grid rules
            [obj] = obj.CollocationGrid();           % Collocation grid
            [obj] = obj.QuadWeights();               % Quadrature weigths
            [obj] = obj.DiffMatrix();                % Differentiation matrix
            obj.J = 1;                               % Jacobian domain transformation 
        end

        % Particular methods
        function [obj] = CollocationGrid(obj)
        end

        function [obj] = DiffMatrix(obj)
            obj.D = [];
        end

        function [obj] = QuadWeights(obj)
        end

        function [t, dt] = Domain(t0, tf, tau)
            t  = (tf-t0) * obj.J * tau; 
            dt = (tf-t0) * obj.J * ones(1,length(tau));
        end
    end   
end