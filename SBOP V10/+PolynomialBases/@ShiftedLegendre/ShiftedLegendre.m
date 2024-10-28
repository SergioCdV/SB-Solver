classdef ShiftedLegendre < PolynomialBases.AbstractBasis

    properties 
    end

    methods         
         [Pn] = basis(obj, order, u);
         [B] = derivative(obj, order, u, degree);
    end
end