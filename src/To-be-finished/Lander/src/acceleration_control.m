%% Project: 
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - vector 
%         - scalar tf, the final time of flight

% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(mu,C,tf,a)
    % Ellipsoidal acceleration
        % Preallocation 
    g = zeros(3,size(C,2)); 

    for i = 1:3
        fun = @(x)( ((a(1)^2+x).*(a(2)^2+x).*(a(3)^2+x)).^(-0.5).*(a(i)^2+x).^(-1) );

        for j = 1:size(C,2)
            % Compute the confocal distance 
            eq = @(x)(C(1,j)^2/(a(1)^2+x)+C(2,j)^2/(a(2)^2+x)+C(3,j)^2/(a(3)^2+x)-1);
            eps = fsolve(eq, a(1));

            % Acceleration 
            L = integral(@(x)fun(x), eps, Inf);
            g(i,j) = -mu*C(i,j)*L;
        end
    end


    % Compute the control vector as a dynamics residual
    u = C(7:9,:)-tf^2*g;
end
