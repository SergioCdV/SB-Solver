classdef Chebyshev < PolynomialBases.AbstractBasis

    properties 
        kind
    end

    methods
        function obj = Chebyshev(myKind)
            if (exist('myKind', 'var'))
                obj.kind = myKind;
            else
                obj.kind = 'first';
            end

        end
         
         [tau, w, D] = weights(N)
         [Pn] = basis(obj, order, u);
         [B] = derivative(kind, order, u, degree);
         [y] = nodes(N);
    end
    
end