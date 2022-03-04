%% Project: 
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar amax, the maximum acceleration allowed for the
%           spacecraft
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - string basis, specifying the polynomial collacation basis

% Outputs: - scalar tfapp, the initial approximation of the time of flight
%          - array Papp, the initial estimation of the boundary control
%            points
%          - array Bapp, the Bernstein polynomials basis in use
%          - array Capp, the initial estimation of the spacecraft state vector

function [Papp, Bapp, Capp] = initial_approximation(mu, tf, tau, n, initial, basis)
    % Approximation order in the Bernstein curve
    n_init = 3; 

    switch (basis)
        case 'Bernstein'
                % Initial estimate of control points (using the boundary conditions)
                Papp = boundary_conditions(mu, tf, n, initial, basis);

                % Bernstein polynomial basis
                Bapp = [bernstein_basis(n_init,tau); bernstein_derivative(n_init,tau,1); bernstein_derivative(n_init,tau,2)];

        case 'Orthogonal Bernstein'
                % Initial estimate of control points (using the non-orthonormal boundary conditions)
                Papp = boundary_conditions(mu, tf, n, initial, basis);

                % Bernstein polynomial basis
                Bapp = [OB_basis(n_init,tau); OB_derivative(n_init,tau,1); OB_derivative(n_init,tau,2)];
        otherwise
            error('No valid collocation polynomial basis has been selected')
    end

    % State vector approximations
    Capp = kron(eye(3),Papp)*Bapp;
end