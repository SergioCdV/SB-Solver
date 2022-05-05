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
        case 'Sundman'
            tau = linspace(0,1,m);
        case 'Normal'
            sigma = 1;                          
            pd = makedist('Normal');
            pd.sigma = sigma;
            xpd = linspace(-3,3,m);
            tau = cdf(pd,xpd);
        case 'Random'
            tau = rand(1,m);
            tau = sort(tau);
        case 'Gauss-Lobatto'
            i = 1:m;
            tau = -cos((i-1)/(m-1)*pi);
            tau = (tau-tau(1))/(tau(end)-tau(1));
        case 'Legendre-Gauss'
            tau = LG_nodes(-1,1,m);
        case 'Chebyshev'
            i = m:-1:1;
            tau = cos((2*i-1)/(2*m)*pi);
        case 'Orthonormal Bezier'
            tau = OB_nodes(0,1,m);
        otherwise
            error('An appropriate time array distribution must be specified')
    end
end