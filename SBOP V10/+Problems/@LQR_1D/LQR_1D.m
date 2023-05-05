%% SBOPT %%
% Date: 05/05/2023

%% 1D LQR %% 
% Implementation of the 1D finite horizon LQR transcription

classdef LQR_1D < Problems.AbstractProblem 
    % Fundamental definition of the problem
    properties  
    end

    methods 
        % Constructor 
        function [obj] = LQR_1D(myInitial, myFinal, myDerDeg, myStateDim, myControlDim)
            obj.initial = myInitial;                % Initial boundary conditions 
            obj.final = myFinal;                    % Final boundary conditions
            obj.DerDeg = myDerDeg;                  % Degree of the derivative              
            obj.StateDim = myStateDim;              % State dimension 
            obj.ControlDim = myControlDim;          % Control dimension 
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
        [LB, UB] = BoundsFunction();
    end
end