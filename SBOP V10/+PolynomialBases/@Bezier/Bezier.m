classdef Bezier < PolynomialBases.AbstractBasis

    properties 
    end

    methods         
         [Pn] = basis(obj, order, u);
         [B] = derivative(obj, order, u, degree);
         [M] = LB_tmatrix(obj, n);
         [M] = DB_tmatrix(obj, n, order);
    end

end