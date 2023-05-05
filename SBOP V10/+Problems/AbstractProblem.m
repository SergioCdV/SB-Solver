%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Dynamics class %% 
% Class implementation of an Optimal Problem function 

classdef AbstractProblem 
    % Fundamental definition of the problem
    properties  
        % State dynamics
        initial;        % Initial boundary conditions 
        final;          % Final boundary conditions
        DerDeg;         % Degree of the derivative              
        StateDim;       % State dimension 
        ControlDim;     % Control dimension

        % General parameters 
        Params;         % General parameters
    end

    methods 
        % Construction 
        function [obj] = AbstractProblem()
            obj = obj.Check();
        end

        % Problem transcription
        [s0, sf] = BoundaryConditions(initial, final, beta, t0, tf);
        [u] = ControlFunction(params, beta, t0, tf, t, s);
        [M, L] = CostFunction(params, beta, t0, tf, s, u);
        [A, b, Aeq, beq] = LinConstraints(beta, P);
        [c, ceq] = NlinConstraints(params, beta, t0, tf, tau, s, u);
        [beta, t0, tf] = InitialGuess(params, initial, final);
    end

    methods (Static)
        % Problem transcription
        [LB, UB] = BoundsFunction();

        % Check object
        function [obj] = Check(obj)
            % Check the dimensionality of the dynamics 
            if (size(obj.initial,1) ~= obj.StateDim * (obj.DerDeg) || size(obj.final,1) ~= obj.StateDim * (obj.DerDeg)) 
                error('Supplied boundary conditions are not of Cauchy type.');
            end

            if (length(obj.PolOrder) ~= obj.StateDim)
                warning('The input polynomial order vector mismatches the state dimension...'); 
                obj.PolOrder = [obj.PolOrder min(obj.PolOrder)*ones(1,obj.StateDim-length(obj.PolOrder))].';
            elseif (size(obj.PolOrder,1) ~= obj.StateDim)
                obj.PolOrder = obj.PolOrder.';
            end

            % Check the order of the dynamics 
            if (obj.DerDeg < 1 || obj.DerDeg > 2)
                error('The problem dynamics are not well-modelled. ODE up to 2 order are supported.');
            end
        end
    end

end