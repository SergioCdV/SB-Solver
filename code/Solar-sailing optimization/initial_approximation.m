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

% Outputs: - scalar tfapp, the initial approximation of the time of flight
%          - array Papp, the initial estimation of the boundary control
%            points
%          - array Bapp, the Bernstein polynomials basis in use
%          - array Capp, the initial estimation of the spacecraft state vector

function [tfapp, Papp, Bapp, Capp] = initial_approximation(mu, r0, amax, tau, initial, final)
    % Approximation order in the Bernstein curve
    n_init = 3; 

    % Initial guess for transfer time (as Hohmann transfer)
    tfapp = 2*abs(sqrt(mu/final(1))*(sqrt(2*initial(10)/(final(10) + initial(10)))*(1-final(10)/initial(10))+sqrt(final(10)/initial(10))-1))/amax;

    % Initial estimate of control points (using the boundary conditions)
    Papp = boundary_conditions(initial, final, tfapp, n_init, r0, 'Bernstein');

    % Bernstein polynomial basis
    Bapp = [bernstein_basis(n_init,tau); bernstein_derivative(n_init,tau,1); bernstein_derivative(n_init,tau,2)];

    % State vector approximations
    n = n_init+1;
    Capp = [Papp*Bapp(1:n,:); Papp*Bapp(n+1:2*n,:); Papp*Bapp(2*n+1:3*n,:)];
end