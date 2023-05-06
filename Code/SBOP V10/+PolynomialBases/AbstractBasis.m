%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Dynamics class %% 
% Class implementation of an Optimal Problem function 

classdef (Abstract) AbstractBasis
    
    methods 
         [Pn] = basis(obj, order, u);
         [B] = derivative(obj, order, u, degree);
    end
end