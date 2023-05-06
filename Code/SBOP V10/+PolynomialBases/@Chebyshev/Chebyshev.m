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
         
         [Pn] = basis(obj, order, u);
         [B] = derivative(obj, order, u, degree);
    end
    
end