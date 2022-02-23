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

function [tfapp, Papp, Bapp, Capp] = initial_approximation(mu, r0, amax, tau, initial, final, basis)
    % Approximation order in the Bernstein curve
    n_init = 3; 

    % Initial guess for transfer time (as Hohmann transfer)
    a = (norm(initial([1 3]))+norm(final([1 3])))/2;
    dE = (1/2)*(norm(final([4 6]))^2 - norm(initial([4 6]))^2)-mu*(1/norm(final([1 3]))-1/norm(initial([1 3])));
    tfapp = abs(dE)/(amax*a*sqrt(mu/a^3));

    switch (basis)
        case 'Bernstein'
                % Initial estimate of control points (using the boundary conditions)
                Papp = boundary_conditions(initial, final, tfapp, n_init, r0, basis);

                % Bernstein polynomial basis
                Bapp = [bernstein_basis(n_init,tau); bernstein_derivative(n_init,tau,1); bernstein_derivative(n_init,tau,2)];
        case 'Orthogonal Bernstein'
                % Linear system of interest 
                b = [initial(1:6)./[r0;1;r0;tfapp/r0;tfapp;tfapp/r0]; final(1:6)./[r0;1;r0;tfapp/r0;tfapp;tfapp/r0]];    % Boundary conditions
                b = [b([1 7 4 10]); b([2 8 5 11]); b([3 9 6 12])];
                A = [OB_basis(n_init,[0 1]).'; OB_derivative(n_init,[0 1],1).'];                                         % Linear system
                A = kron(eye(3),A);

                % Initial estimate of control points (using the boundary conditions)
                Papp = A\b; 
                Papp = reshape(Papp, [n_init+1 3]).';

                % Bernstein polynomial basis
                Bapp = [OB_basis(n_init,tau); OB_derivative(n_init,tau,1); OB_derivative(n_init,tau,2)];
        otherwise
            error('No valid collocation polynomial basis has been selected')
    end

    % State vector approximations
    Capp = kron(eye(3),Papp)*Bapp;
end