%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Constants 
    mu = params(1);             % Gravitational parameter
    COE = params(2:7);          % Orbit elements
    T = params(17:19);          % Control authority
    omega0 = params(21:23);     % Initial angular velocity
    omegaf = params(24:26);     % Final angular velocity
    omega_max = params(20);     % Maximum angular velocity

    n_ECI = params(27:29);      % Pointing direction
    ng = params(30:32);         % Body frame vector to which we shall point
    thetanp = params(33);       % Geo-pointing tolerance
    
    % Anomaly space 
    v = (tf-t0) * 0.5 * (1+tau) + COE(end);

    % Transformation to the ECI frame 
    E = OrbitalDynamics.euler_matrix(COE);      % Euler matrix
    D = [0 1 0; 0 0 -1; -1 0 0];                % Alignment matrix

    % Magnetic field across the time span in the LVLH frame
    B = OrbitalDynamics.dipole_model(mu, COE, v);

    for i = 1:size(B,2)
        % Transformation to the perifocal frame 
        R = [cos(v(i)) sin(v(i)) 0; -sin(v(i)) cos(v(i)) 0; 0 0 1];
        
        % Transformation to the body frame
        S = [0 -s(3,i) s(2,i); s(3,i) 0 -s(1,i); -s(2,i) s(1,i) 0];
        A = (s(4,i)^2-s(1:3,i).'*s(1:3,i)) * eye(3) + 2 * s(1:3,i) * s(1:3,i).' - 2 * s(4,i) * S;
        B(:,i) = (A * (D * R * E).') * B(:,i);
    end

    ub = u; 
    for i = 1:length(tau)
        ub(:,i) = cross(B(:,i), u(:,i)) / norm(B(:,i))^2;
    end

    % Angular velocity 
    omega = 2 * [s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
                -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
                 s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];

    % Inequality constraints
    sat = abs(ub) - repmat(T, 1, size(ub,2));

    c = [reshape(sat, 1, []) ...                                    % Torque saturation
         dot(omega, omega, 1)-omega_max^2 ...                       % Angular velocity tolerance
         dot(s(1:4,:), s(1:4,:),1)-ones(1,size(u,2)) ...            % Quaternion norm constraint
         acos(abs(dot(n_ECI, A.' * ng)))-thetanp];                  % Geo-pointing requirement

    % Fitting of the control law
    for i = 1:length(tau)
        ub(:,i) = cross(ub(:,i), B(:,i));
    end

    c = [c dot(ub-u, ub-u, 1) - 1e-8];

    % Equality constraints
    ceq = [omega(1:3,1) - omega0; ...
           omega(1:3,end)- omegaf].';

%            cos( tf-t0 )-cos(COE(end)); ...
%            sin( tf-t0 )-sin(COE(end))
end