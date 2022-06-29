%% Project: 
% Date: 16/04/22

%% Collocation grid %%
% Function to compute the collocation grid

% Inputs: - scalar m, the number of collocation
%         - string method, the type of grid to use
%         - string mode, to select if the sampling grid is based on the
%           the intersection of low-order nodes or the nodes of a
%           high-order polynomial

% Outputs: - vector tau, the collocation points to be used 

function [tau] = collocation_grid(m, method, mode)
    switch (mode)
        case 'Intersection'
            tau = zeros(1,sum(2:m));
            for i = 2:m
                tau(1+sum(2:i-1):sum(2:i)) = sampling_grid(i,method);
            end
            tau = unique(tau);
            tau = sort(tau);
        otherwise
            tau = sampling_grid(m,method);
    end
end

%% Auxiliary functions 
% Sampling grid computation
function [tau] = sampling_grid(m,method)
    switch (method)
        case 'Linear'
            tau = linspace(0,1,m);
        case 'Normal'
            sigma = 1;                          
            pd = makedist('Normal');
            pd.sigma = sigma;
            xpd = linspace(-sigma,sigma,m-2);
            tau = cdf(pd,xpd);
            tau = [0 tau 1];  
        case 'Random'
            tau = rand(1,m);
            tau = sort(tau);
        case 'Legendre'
            tau = LG_nodes(m);
        case 'Chebyshev'
            tau = CH_nodes(m);
        case 'Laguerre'
            tau = LR_nodes(m,0);
        case 'Hermite'
            % Depricated
            tau = HT_nodes(m);
        case 'Orthogonal Bernstein'
            % Depricated
            tau = OB_nodes(m);
        case 'Regularized'
            tau = linspace(0,1,m);
        otherwise
            error('An appropriate time array distribution must be specified')
    end
    if (size(tau,1) ~= 1)
        tau = tau.';
    end
end