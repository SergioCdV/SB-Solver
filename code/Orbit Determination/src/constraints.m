%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - array measurements, an mx4 matrix of measurements in the form
%           of epoch | 3D vector measurement
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string cost_function, the minimization policy to follow
%         - string basis, the polynomial basis to be used
%         - string method, the parameter distribution to be used

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, T, initial, final, measurements, n, x, B, cost_function, basis, method)
    % Extract the optimization variables
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % Final time of flight 
    N = floor(x(end));                                  % Optimal number of revolutions

    tol = 1;                                            % Measurements error tolerance

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Radius constraints
    R = sqrt(C(1,:).^2+C(3,:).^2);
    r0 = sqrt(initial(1)^2+initial(3)^2);

    if (~isempty(final))
        rf = sqrt(final(1)^2+final(3)^2);
    else
        rf = r0;
    end

    % Equalities 
    ceq = [];

    % Inequality (control authority)
    switch (cost_function)
        case 'Dynamics residual'
            epochs = measurements(1,:)/tf; 
            switch (method)
                case 'Regularized'
                    options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
                    [~, epochs] = ode45(@(t,s)Sundman_transformation(basis,n,P,t,s), epochs, 0, options);
                    epochs = epochs.';
                case 'Chebyshev'
                    epochs = 2*epochs-1;
                case 'Legendre'
                    epochs = 2*epochs-1;
                case 'Laguerre'
                    epochs = 2*epochs-1;
            end
    
            % State evolution
            B = state_basis(n, epochs, basis);
            C = evaluate_state(P,B,n);
        
            % Compute the residuals 
            M = cylindrical2cartesian(C(1:3,:),true);
            M(1:3,:) = M(1:3,:)./sqrt(M(1,:).^2+M(2,:).^2+M(3,:).^2);
            e = measurements(2:4,:)-M(1:3,:);
            r = dot(e,e,1);
            r = mean(r)+3*std(r);

            c = [R-2*max([r0 rf]) r-tol];
        otherwise
            % Control input 
            u = acceleration_control(mu,C,tf,method);

            switch (method)
                case 'Regularized'
                    c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-r.^2.*tf^2*T*ones(1,size(u,2)) R-2*max([r0 rf])]; 
                otherwise
                    c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-tf^2*T*ones(1,size(u,2)) R-2*max([r0 rf])];
            end
    end
end

%% Auxiliary function
% Compute the Sundman transformation 
function [ds] = Sundman_transformation(basis, n, P, t, s)
    B = state_basis(n, s, basis);
    C = evaluate_state(P,B,n);
    ds = (C(1,:).^2+C(3,:).^2).^(-1/2);
end