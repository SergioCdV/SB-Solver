classdef (Abstract) AbstractBasis
    
    methods 
         [tau, w, D] = weights(N);
         [Pn] = basis(order, u);
         [B] = derivative(order, u, degree);
         [y] = nodes(N);
    end
end