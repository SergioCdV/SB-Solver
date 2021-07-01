function [tfapp, Papp, Bapp, Capp] = initial_approximation(initial, final, tau, m, mu, amax, r0)

n_init = 3; % order of bezier functions

%% Initial guess for transfer time (hohmann transfer
tfapp = 2*abs(sqrt(mu/final.pos(1))*(sqrt(2*initial.a/(final.a + initial.a))*(1-final.a/initial.a)+sqrt(final.a/initial.a)-1))/amax;


%% Initial estimate of control points
Papp = boundary_conditions(initial, final, tfapp, n_init, r0);


%% initial estimate of bernstein polynomials
Bapp = zeros(3,m,n_init+1); % (coordinates, discretisation, and control point)

for i = 0:2 % cicle through 0th, 1st, 2nd derivatives
    for j=0:n_init % for each of the control points
       Bapp(i+1,:,j+1) = bernstein(n_init,j,tau,i);
    end
end



%% Initial approximation data values (coordinate evolution)

% first index indicates the derivative order (0th, 1st, 2nd)
% second index indicates the discretization step
% third index indicates the coordinate (rho, theta, z)
Capp = zeros(3,m,3);

for i=1:3 % for coordinates
    for j=1:3 % for derivatives
        Capp(j,:,i) = squeeze(Bapp(j,:,:))*Papp(:,i);
    end
end

end