%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - dynamics, string specifying the independent variable
%           determining the dynamics of the problem
%         - vector tau, the collocation points to be used 
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collacation basis

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - initial estimation of the control law
%          - scalar tfapp, the initial initial time of flight

function [Papp, Uapp, tfapp] = initial_approximation(tau, tfapp, initial, final, basis)
    % Constants 
    n_init = [3 3 3];

    % Preliminary number of revolutions
    dtheta = final(2)-initial(2);
    if (dtheta < 0)
        dtheta = dtheta + 2*pi; 
    end
    
    Napp = ceil( (dtheta+tfapp*0.5*(initial(4)+final(4)) ) / (2*pi) );
    if (Napp <= 0)
        Napp = 1;
    end 
    
    % New initial TOF
    tfapp = tfapp*Napp;

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(length(initial)/2, max(n_init)+1);
    [B, ~] = state_basis(n_init, tau, basis);
    Uapp = evaluate_state(Papp, B, n_init);
end