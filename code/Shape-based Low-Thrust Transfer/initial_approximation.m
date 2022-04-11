%% Project: 
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - vector tau, the collocation points to be used 
%         - vector n, the degree of approximation of each state variable 
%         - scalar T, the maximum acceleration allowed for the
%           spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collacation basis

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - array Bapp, the Bernstein polynomials basis in use
%          - array Capp, the initial estimation of the spacecraft state vector
%          - scalar tfapp, the initial approximation of the time of flight
%          - scalar N, the estimated number of revolutions needed

function [Papp, Bapp, Capp, tfapp, N] = initial_approximation(mu, tau, n, T, initial, final, basis)
    % Approximation order in the Bernstein curve
    n_init = 3; 

    % Approximation of the time of flight 
    r0 = norm(initial([1 3]));              % Initial radius
    v0 = norm(initial(4:6));                % Initial velocity
    rf = norm(final([1 3]));                % Final radius
    vf = norm(final(4:6));                  % Final velocity
    a = (r0+rf)/2;                          % Transfer semimajor axis
    dE = (vf^2-v0^2)/2-mu*(1/rf-1/r0);      % Change of energy
    tfapp = abs(dE)/(T*sqrt(mu/a));         % Time of flight
    tfapp = 2*pi*(800/365);

    % Preliminary number of revolutions 
    N = floor(initial(2)-final(2)+tfapp*(initial(4)+final(4))/(2*pi))-1;
    if (N < 0)
        N = 0;
    end

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Papp = boundary_conditions(tfapp, n, initial, final, N, basis);

    % Generate the polynomial basis
    switch (basis)
        case 'Bernstein'
                % Bernstein polynomial basis
                Bapp = [bernstein_basis(n_init,tau); bernstein_derivative(n_init,tau,1); bernstein_derivative(n_init,tau,2)];

        case 'Orthogonal Bernstein'
                % Bernstein polynomial basis
                Bapp = [OB_basis(n_init,tau); OB_derivative(n_init,tau,1); OB_derivative(n_init,tau,2)];
        otherwise
            error('No valid collocation polynomial basis has been selected')
    end

    % State vector approximations
    Capp = kron(eye(size(Bapp,1)/size(Papp,2)),Papp)*Bapp;
end