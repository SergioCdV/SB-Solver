%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    g = params(1) * s(1:3,:); 
    options = optimset('Display', 'off');

    % Gravity field
    for i = 1:size(s,2)
        % Compute the confocal distance 
        xi_fun = @(x)(s(1,i)^2/(params(2)^2+x^2) + s(2,i)^2/(params(3)^2+x^2) + s(3,i)^2/(params(4)^2+x^2) - 1);
        xi = fsolve(xi_fun, params(1), options);

        % Integrate the work
        alpha = zeros(3,1);
        for j = 1:3
            fun = @(alpha)( 1./(sqrt((params(2)^2 + alpha ).*(params(3)^2 + alpha).*(params(4)^2 + alpha)) .* (params(1+j)^2 + alpha)) );
            alpha(j) = integral(fun, xi, Inf);
        end

        g(:,i) = g(:,i) .* alpha;
    end
  
    % Compute the control vector as a dynamics residual
    u = s(7:9,:) - g;
end