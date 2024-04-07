%% SBOPT %%
% Date: 05/05/2023

%% Robot problem %% 
% Implementation of the 1D robot

classdef RobotKinematics < Problems.AbstractProblem 
    % Fundamental definition of the problem
    properties  
    end

    methods 
        % Constructor 
        function [obj] = RobotKinematics(myInitial, myFinal, myDerDeg, myStateDim, myControlDim, myParams)
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
        [s0, sf] = BoundaryConditions(obj, initial, final, beta, t0, tf);
        [u] = ControlFunction(obj, params, beta, t0, tf, t, s);
        [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u);
        [A, b, Aeq, beq] = LinConstraints(obj, params, beta, P);
        [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u);
        [beta, t0, tf] = InitialGuess(obj, params, initial, final);
        [LB, UB] = BoundsFunction(obj);

        function [obj] = Check(obj)
            obj = Check@Problems.AbstractProblem(obj);
        end
    end

    methods (Static)
        [T, J] = Kinematics(obj, n, transformation, s);
        [T, R] = ur3_dkinematics(obj, i, q);
        [dm] = collision_constraint(r1, A1, v1, r2, A2, v2, iterations);
        [v] = body_vertex(index);
    end
end