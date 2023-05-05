classdef Legendre < PolynomialBases.AbstractBasis

    properties 
    end

    methods         
         [tau, w, D] = weights(N)
         [Pn] = basis(obj, order, u);
         [B] = derivative(order, u, degree);
         [y] = nodes(N);
    end
    
end