classdef OBezier < PolynomialBases.AbstractBasis

    properties 
        kind
    end

    methods         
         [tau, w, D] = weights(N)
         [Pn] = basis(obj, order, u);
         [B] = derivative(kind, order, u, degree);
         [y] = nodes(N);
    end

end