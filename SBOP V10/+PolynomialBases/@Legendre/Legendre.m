classdef Legendre < PolynomialBases.AbstractBasis

    properties 
    end

    methods         
         [Pn] = basis(obj, order, u);
         [B] = derivative(obj, order, u, degree);
         [C] = modal_projection(obj, S);
    end
end