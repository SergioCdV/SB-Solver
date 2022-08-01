%% Project: 
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - vector tau, the collocation points to be used 
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collacation basis

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - array Capp, the initial estimation of the spacecraft state vector

function [Papp, Capp] = initial_approximation(tau, tfapp, initial, final, basis)
    % Generate the polynomial basis
    n_init = 3; 
    n = [n_init n_init n_init];
    Bapp = state_basis(n, tau, basis);

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(length(initial)/2, max(n_init)+1);  
    Papp = boundary_conditions(tfapp, n, initial, final, Papp, Bapp, basis);

    % State vector approximations
    Capp = evaluate_state(Papp, Bapp, n);
end