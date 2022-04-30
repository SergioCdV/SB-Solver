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
%          - scalar Napp, the estimated number of revolutions needed
%          - scalar tfapp, the initial initial time of flight

function [Papp, Capp, Napp, tfapp] = initial_approximation(tau, tfapp, initial, final, basis)
    % Approximation order in the Bernstein curve
    n_init = 3; 

    % Preliminary number of revolutions 
    Napp = floor(initial(2)-final(2)+tfapp*(initial(4)+final(4))/(2*pi))-1;
    if (Napp <= 0)
        Napp = 1;
    end

    % New initial TOF
    tfapp = tfapp*Napp;

    % Generate the polynomial basis
    n = [n_init n_init n_init];
    Bapp = state_basis(n,tau,basis);

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = boundary_conditions(tfapp, n_init, initial, final, Napp, [], Bapp, basis);

    % State vector approximations
    Capp = evaluate_state(Papp, Bapp, n);
end