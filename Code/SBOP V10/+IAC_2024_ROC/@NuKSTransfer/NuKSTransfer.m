%% SBOPT %%
% Date: 16/08/2023

%% KS orbital transfer problem %% 
% Implementation of optimal low-thrust 3D orbital transfers in KS
% coodinates and true anomaly regularization

classdef NuKSTransfer < Problems.AbstractProblem 
    % Fundamental definition of the problem
    properties  
    end

    methods 
        % Constructor 
        function [obj] = NuKSTransfer(myInitial, myFinal, myDerDeg, myStateDim, myControlDim, myParams)
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
        % KS dynamics
        [L] = KS_matrix(u);
        [L] = IKS_matrix(u);
        [l] = bilinear_function(u,v);
        [u] = u_mapping(r);
        [S] = state_mapping(x, direction);
    end
end