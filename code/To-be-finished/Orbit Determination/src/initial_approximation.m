%% Project: 
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - sampling_distribution, string specifying the sampling grid
%           distribution
%         - vector tau, the collocation points to be used 
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collacation basis
%         - string dynamics, the parametrization of the dynamics
%           vectorfield to be used

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - array Capp, the initial estimation of the spacecraft state vector
%          - scalar Napp, the estimated number of revolutions needed
%          - scalar tfapp, the initial initial time of flight

function [Papp, Capp, Napp, tfapp] = initial_approximation(dynamics, tau, tfapp, initial, final, basis)
    % Generate the polynomial basis
    n_init = 3; 
    n = repmat(n_init,1,3);
    Bapp = state_basis(n,tau,basis);

    if (~isempty(final))
        % Preliminary number of revolutions
        dtheta = final(2)-initial(2);
        if (dtheta < 0)
            dtheta = dtheta + 2*pi; 
        end
        
        Napp = ceil( (dtheta+tfapp*0.5*(initial(4)+final(4)) ) / (2*pi) );
        if (Napp <= 0)
            Napp = 1;
        end
    else
        Napp = 0; 
    end 
        
    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = zeros(length(initial)/2, n_init+1);  
    Papp = boundary_conditions(tfapp, n, initial, final, Napp, Papp, Bapp, basis);

    % State vector approximations
    Capp = evaluate_state(Papp, Bapp, n);

    % Time-regularized solution 
    switch (dynamics)
        case 'Regularized'
            % Arc-length regularization
            r = sqrt(Capp(1,:).^2+Capp(3,:).^2);
            tfapp = tfapp*trapz(tau, r.^(-1));  
            Papp = boundary_conditions(tfapp, n, initial, final, Napp, Papp, Bapp, basis);
        
            % State vector approximations
            Capp = evaluate_state(Papp, Bapp, n);
        otherwise
    end
end