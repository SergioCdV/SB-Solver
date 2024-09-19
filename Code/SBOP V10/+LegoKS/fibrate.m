%% Fibrate %%
% Function to compute the a fiber associated to a given KS vector

% Inputs: - vector u, one u vector 
%         - vector theta, a distribution of the phase angle

% Outputs: - F, the associated fiber 

function [F] = fibrate(u, theta)
    % Pre-allocation 
    F = zeros( size(u) );

    N = size(u,1) / 4;

    for i = 1:size(theta,2)
        R = LegoKS.GroupAction( theta(i) );
        R = repmat(R, 1, N);
        R = mat2cell( R, 4, 4 * ones(1,N) );
        Rt = blkdiag( R{:} );
        F(:,i) = Rt * u(:,1);
    end
end