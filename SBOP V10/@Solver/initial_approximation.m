%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - object Problem, defining the transcription of the problem of interest
%         - cell array B, containing the polynomial support of the state
%           vector
%         - object Grid, the collocation points object to be used 

% Outputs: - vector beta, the initial estimation of the optimization extra variables
%          - scalar t0, the initial initial time 
%          - scalar tf, the initial estimation of the final time of flight
%          - array P, the initial estimation of the control points
%          - array C, the initial estimation of the spacecraft state vector

function [beta, t0, tf, P, C] = initial_approximation(obj, Problem, B, Grid)
    % Initial guess 
    [beta, t0, tf] = Problem.InitialGuess(Problem.Params, Problem.initial, Problem.final);

    % Initial estimate of state points
    tau = Grid.Domain(t0, tf, Grid.tau);
    P = zeros(Problem.StateDim, max(n_init)+1);  
    P = obj.boundary_conditions(Problem, beta, t0, tf, tau, B, P);

    % State vector approximation as a function of time
    C = obj.evaluate_state(n, L, P, B);
end