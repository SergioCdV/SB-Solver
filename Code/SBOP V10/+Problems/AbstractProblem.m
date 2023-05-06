%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Dynamics class %% 
% Class implementation of an Optimal Problem function 

classdef (Abstract) AbstractProblem 
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
        function [obj] = AbstractProblem(myInitial, myFinal, myDerDeg, myStateDim, myControlDim, myParams)
            obj.initial = myInitial;                % Initial boundary conditions 
            obj.final = myFinal;                    % Final boundary conditions
            obj.DerDeg = myDerDeg;                  % Degree of the derivative              
            obj.StateDim = myStateDim;              % State dimension 
            obj.ControlDim = myControlDim;          % Control dimension
            obj.Params = myParams;                  % Problem parameters
        end

        % Problem transcription
        [s0, sf] = BoundaryConditions(obj, initial, final, beta, t0, tf);
        [u] = ControlFunction(obj, params, beta, t0, tf, t, s);
        [M, L] = CostFunction(obj, params, beta, t0, tf, s, u);
        [A, b, Aeq, beq] = LinConstraints(obj, beta, P);
        [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u);
        [beta, t0, tf] = InitialGuess(obj, params, initial, final);
        [LB, UB] = BoundsFunction(obj);
        
        function [obj] = Check(obj)
            % Check the dimensionality of the dynamics 
            if (size(obj.initial,1) ~= obj.StateDim * (obj.DerDeg) || size(obj.final,1) ~= obj.StateDim * (obj.DerDeg)) 
                error('Supplied boundary conditions are not of Cauchy type.');
            end

            % Check the order of the dynamics 
            if (obj.DerDeg < 1 || obj.DerDeg > 2)
                error('The problem dynamics are not well-modelled. ODE up to 2 order are supported.');
            end
        end
    end
end