%% Project: 
% Date: 16/04/22

%% Collocation grid %%
% Function to compute the collocation grid

% Inputs: - scalar m, the number of collocation
%         - string method, the type of grid to use

% Outputs: - vector tau, the collocation points to be used 

function [tau] = collocation_grid(m, method)
    switch (method)
        case 'Linear'
            tau = linspace(0,1,m);
        case 'Normal'
            sigma = 1;                          
            pd = makedist('Normal');
            pd.sigma = sigma;
            xpd = linspace(-sigma,sigma,m);
            tau = cdf(pd,xpd);
        case 'Random'
            tau = rand(1,m);
            tau = sort(tau);
        case 'Legendre'
            tau = LG_nodes(m);
        case 'Chebyshev'
            tau = CH_nodes(m);
        case 'Laguerre'
            tau = LR_nodes(0,1,m);
        case 'Orthogonal Bernstein'
            tau = OB_nodes(m);
        case 'Sundman'
            tau = linspace(0,1,m);
        otherwise
            error('An appropriate time array distribution must be specified')
    end
end