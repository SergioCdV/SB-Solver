%% SBOPT %%
% Date: 05/05/2023

%% 1D Minimum time trajectory %% 
% Implementation of the 1D minimum time trajectory transcription

classdef MinTime_1D < Problems.AbstractProblem 
    % Fundamental definition of the problem
    properties  
    end

    methods 
        % Constructor 
        function [obj] = MinTime_1D(myInitial, myFinal, myDerDeg, myStateDim, myControlDim, myParams)
            super_arguments{1} = myInitial;
            super_arguments{2} = myFinal;
            super_arguments{3} = myDerDeg;
            super_arguments{4} = myStateDim;
            super_arguments{5} = myControlDim;

            if (exist('myParams', 'var'))
                super_arguments{6} = myParams;
            else
                super_arguments{6} = [];
            end

            obj@Problems.AbstractProblem(super_arguments{:});

            % Check the problem definition
            obj = obj.Check();
        end

        % Problem transcription
        [s0, sf] = BoundaryConditions(initial, final, beta, t0, tf);
        [u] = ControlFunction(params, beta, t0, tf, t, s);
        [M, L] = CostFunction(params, beta, t0, tf, s, u);
        [A, b, Aeq, beq] = LinConstraints(beta, P);
        [c, ceq] = NlinConstraints(params, beta, t0, tf, tau, s, u);
        [beta, t0, tf] = InitialGuess(params, initial, final);

        function [obj] = Check(obj)
            obj = Check@Problems.AbstractProblem(obj);
        end
    end

    methods (Static)
        [LB, UB] = BoundsFunction();
    end
end