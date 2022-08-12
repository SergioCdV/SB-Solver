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
%          - array Capp, the initial estimation of the spacecraft state vector
%          - scalar Napp, the estimated number of revolutions needed
%          - scalar tfapp, the initial initial time of flight

function [Papp, Capp, Napp, tfapp] = initial_approximation(dynamics, tau, tfapp, initial, final, basis)
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

    % Generate the polynomial basis
    n_init = repmat(3, [1 3]);
    Bapp = state_basis(n_init, tau, basis);

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(length(initial)/2, max(n_init)+1);  
    Papp = boundary_conditions(tfapp, n_init, initial, final, Napp, Papp, Bapp, basis);

    % State vector approximation as a function of time
    Capp = evaluate_state(Papp, Bapp, n_init);

    % Time-regularized solution 
    switch (dynamics)
        case 'Sundman'
            % Arc-length regularization
            r = sqrt(Capp(1,:).^2+Capp(3,:).^2);
            tfapp = trapz(tau, tfapp./r);  
            Papp = zeros(length(initial)/2, max(n_init)+1);  
            Papp = boundary_conditions(tfapp, n_init, initial, final, Napp, Papp, Bapp, basis);
        
            % State vector approximation as a function of s
            Capp = evaluate_state(Papp, Bapp, n_init);
            
        otherwise
    end
end